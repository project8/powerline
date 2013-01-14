#include "correline_utils.hh"
#include <math.h>


//-------------Background------------

bool Background::load(string fname)
{
	ifstream fin(fname.c_str());
	if(!fin.good()) {
		cerr << "unable to open file: " << fname << endl;
		return false;
	}
	//count bg size
	bg_size=0;
	string line;
	while(getline(fin,line))
		if(line[0]!='#')
			bg_size++;
	background_avg=new fftwf_complex[bg_size];
	background_stdev=new fftwf_complex[bg_size];
	background_stdev_invert=new fftwf_complex[bg_size];
	fin.clear();
	fin.seekg(ios_base::beg);
	int on_bg=0;
	double scale=1;
	double f;
	while(getline(fin,line)) {
		if(line[0]!='#')
		{
			stringstream ss(line);
			ss >> f;
			ss >> background_avg[on_bg][0];
			ss >> background_avg[on_bg][1];
			ss >> background_stdev[on_bg][0];
			ss >> background_stdev[on_bg][1];
			on_bg++;
		} else {
//			cerr << "line: |" << line << "|" << endl;
			if(line.compare(0,7,"#scale:")==0) {
				stringstream ss(line);
				string blah;
				ss >> blah >> scale;
//				cerr << "scale is " << scale << endl;
			}
		}
	}
	fin.close();
	for(int i=0;i<bg_size;i++) {
		for(int j=0;j<2;j++)
			if(background_stdev[i][j]!=0)
				background_stdev_invert[i][j]=1/background_stdev[i][j];
			else
				background_stdev_invert[i][j]=1;
	}

	//remove normalization 
	for(int i=0;i<bg_size;i++)
		for(int j=0;j<2;j++) {
			background_avg[i][j]/=scale;
			background_stdev[i][j]/=scale;
			background_stdev_invert[i][j]*=scale;
		}
	return true;
}

Background::~Background()
{
	delete background_avg;
	delete background_stdev;
	delete background_stdev_invert;
}

//------------------ConvolutionMap-------------
int ConvolutionMap::get_t_extent()
{
	int tmax=0;
	for(int i=0;i<conv_length;i++)
		if(conv_ts[i]>tmax)
			tmax=conv_ts[i];
	return tmax;
}

int ConvolutionMap::get_f_extent()
{
	int fmax=0;
	for(int i=0;i<conv_length;i++)
		if(conv_fs[i]>fmax)
			fmax=conv_fs[i];
	return fmax;
}

bool ConvolutionMap::load(string fname) 
{
	ifstream fin(fname.c_str());
	if(!fin.good()) {
		cerr << "unable to open file: " << fname << endl;
		return false;
	}
	conv_length=0;
	string line;
	while(getline(fin,line)) {
		if(line[0]!='#')
		{
			stringstream ss(line);
			ss >> conv_ts[conv_length] >> conv_fs[conv_length];
			conv_length++;
		}
	}
	fin.close();
	return true;
}

void ConvolutionMap::convolve(Waterfall *input,Waterfall *output)
{
	float normalization=1/sqrt((float)conv_length);
	int tl=get_t_extent();
	int fl=get_f_extent();
	for(int i=0;i<input->npoints_t-tl;i++) {
		for(int j=0;j<input->npoints_f-fl;j++) {
			int oindex=output->getIndex(j,i);
			output->data[oindex][0]=0;
			output->data[oindex][1]=0;
			for(int k=0;k<conv_length;k++) {
				int index=output->getIndex(j+conv_fs[k],i+conv_ts[k]);
				output->data[oindex][0]+=input->data[index][0];
				output->data[oindex][1]+=input->data[index][1];
			}
			output->data[oindex][0]*=normalization;
			output->data[oindex][1]*=normalization;
		}
	}
}


//------------------Candidate------------

vector<Candidate> loadCandidateFile(string fname)
{
    long long npoints;
    return loadCandidateFile(fname,npoints);
}

vector<Candidate> loadCandidateFile(string fname,long long &npoints)
{
	vector<Candidate> ret;
	ifstream fin(fname.c_str());
	if(!fin.good()) {
		cerr << "error opening candidate file " << fname << endl;
	}
	string line;
	while(getline(fin,line)) {
		stringstream ls(line);
		if(line[0]=='#') {
		    if(line.compare(0,8,"#npoints")==0) {
			string dummy;
			ls >> dummy >> npoints;
		    }
		    continue;
		}
		Candidate toadd;
		ls >> toadd.frequency >> toadd.time >> toadd.magnitude >> toadd.probability;
		ret.push_back(toadd);
	}
	fin.close();
	return ret;
}

int Candidate::getEventNo(int sampling_rate_mhz,int event_size)
{
	double event_length=((double)event_size)/(sampling_rate_mhz*1e6);
	double event_frac=time/event_length;
	return (int)(floor(event_frac));
}
	

	
bool SortCandidatesByPower::operator()(const Candidate &a,const Candidate &b)
{
	if(a.magnitude<b.magnitude) return true;
	if(a.magnitude>b.magnitude) return false;
	if(a.time<b.time) return true;
	if(a.time>b.time) return false;
	return a.frequency<b.frequency;
}


void PowerCandidateList::push(const Candidate &c)
{
	insert(c);
	if(size()>max_length) erase(begin());
}
	
void PowerCandidateList::saveToFile(string fname,long long npoints)
{
	ofstream fout(fname.c_str());
	fout << "#npoints " << npoints << endl;
	for(PowerCandidateList::iterator it=begin();it!=end();it++) {
		fout << (*it).frequency << " " << (*it).time << " " << (*it).magnitude << " " <<  (*it).probability << endl;
	}
	fout.close();
}
	
PowerCandidateList PowerCandidateList::condenseByDistance(double freqscale,double timescale)
{
	PowerCandidateList ret;
	for(iterator it=begin();it!=end();it++) {
		iterator neighbor=ret.end();
		for(iterator jt=ret.begin();jt!=ret.end();jt++) {
			if(fabs((*jt).frequency-(*it).frequency)<freqscale)
			if(fabs((*jt).time-(*it).time)<timescale)
				neighbor=jt;
		}
		if(neighbor!=ret.end()) {
			if((*neighbor).magnitude<(*it).magnitude) {
				ret.erase(neighbor);
				ret.insert((*it));
			}
		} else
			ret.insert(*it);
	}
	return ret;
}

#include "Waterfall.hh"
#include "correline_utils.hh"
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

string input_eggname;
string background_filename;
string convolution_filename;
string candidate_list_filename;
string prefix;
int tbins=20;
int fbins=20;

Background background;

void print_usage(); //print out instructions
int handle_options(int argc,char *argv[]); //command line options

int main(int argc,char *argv[])
{
	int onindex;
    if(((onindex=handle_options(argc,argv))==-1)||(argc-onindex<0))
	{
		print_usage();
		return -1;
	}
	//----test mandatory inputs provided----
	if(input_eggname=="") {
		cerr << "no egg file name given" << endl;
		return -1;
	}
	if(background_filename=="") {
		cerr << "no background file name given" << endl;
		return -1;
	}
	if(convolution_filename=="") {
		cerr << "no convolution file name given" << endl;
		return -1;
	}
	if(candidate_list_filename=="") {
		cerr << "no candidate list given " << endl;
		return -1;
	}
	//----load files ----
	if(!background.load(background_filename)) {
		return -1;
	}
	vector<Candidate> candidates=loadCandidateFile(candidate_list_filename);
	if(candidates.size()==0) {
		cerr << "no candidates to view" << endl;
		return -1;
	}
	int fft_size=(background.bg_size-1)*2;
    const Monarch *egg=Monarch::OpenForReading(std::string(input_eggname));
	Correlator correlator;
	correlator.init(egg,fft_size);
	//----load the convolution----
	ConvolutionMap cmap;
	cmap.load(convolution_filename);

	//get candidate events
	for(unsigned int i=0;i<candidates.size();i++) {
		candidates[i].event_no=candidates[i].getEventNo(correlator.sampling_rate_mhz,correlator.record_size);
	}
	
	int on_event=0;
	//Start reading egg
	double total_time=0;
    while(egg->ReadRecord()) {
		//see if there are any candidates in this event
		bool candidate_here=false;
		for(unsigned int i=0;i<candidates.size();i++) {
			if(candidates[i].event_no==on_event)
				candidate_here=true;
		}
		if(candidate_here) {
			correlator.process_event(egg);
			for(unsigned int i=0;i<candidates.size();i++) {
				if(candidates[i].event_no==on_event) {
					cout << "on frequency " << candidates[i].frequency << endl;
					cout << "fstep is " << correlator.output_waterfall.freq_step << endl;

					int fbin=correlator.output_waterfall.getFrequencyBin(candidates[i].frequency*1e6);
					int tbin=correlator.output_waterfall.getTimeBin(candidates[i].time-total_time);
					cout << "fbin " << fbin << " tbin " << tbin << endl;
					stringstream onames;
					onames << prefix << "candidate_" << i << ".txt";
	//				/*
					string fname=onames.str();
					ofstream fout(fname.c_str());
					for(int k=fbin-fbins/2;k<fbin+fbins/2;k++) {
					for(int m=tbin-tbins/2;m<tbin+tbins/2;m++) {
					//	int index=correlator.output_waterfall.getIndex(k,m);
					//	double x=correlator.fft_output_1[index][0];
					//	double y=correlator.fft_output_1[index][1];
					//	fout << k << " " << m << " " << x*x+y*y << endl;
						double v=0;
						for(int n=0;n<cmap.conv_length;n++) {
							int nindex=correlator.output_waterfall.getIndex(k+cmap.conv_fs[n],m+cmap.conv_ts[n]);
							double x=correlator.getPowerChannel1(nindex)*background.background_stdev_invert[k+cmap.conv_fs[n]][0];
							double y=correlator.getPowerChannel2(nindex)*background.background_stdev_invert[k+cmap.conv_fs[n]][1];
//							double x=correlator.fft_output_1[nindex][0];
//							double y=correlator.fft_output_1[nindex][1];
							//v+=x*x+y*y;
							v+=0.5*(x+y);
						}
						fout << k << " " << m << " ";
						fout << ((double)(k))*correlator.output_waterfall.freq_step << " " << ((double)(m))*correlator.output_waterfall.time_step+total_time << " ";
						fout << v << endl;
					}
					fout << endl;
					}
					fout.close();
	//				*/
	//				correlator.output_waterfall.saveSubWaterfall(fbin-fbins/2,fbin+fbins/2,tbin-tbins/2,tbin+tbins/2,total_time,onames.str());
				}
			}
		}
		on_event++;
		total_time+=correlator.getEventDuration();
	}
}

void print_usage()
{
	cout << "view_candidate" << endl;
	cout << "prints out a window around a candidate" << endl;
	cout << "options: " << endl;
	cout << "    -i (filename) sets the input egg file MANDATORY" << endl;
	cout << "    -b (filename) sets the ascii background correlation spectrum MANDATORY" << endl;
	cout << "    -l (filename) sets the list of candidates MANDATORY" << endl;
	cout << "    -c (filename) sets the convolution shape" << endl;
	cout << "    -f (integer) sets the number of frequency bins around the candidate to plot" << endl;
	cout << "    -t (integer) sets the number of time bins around the candidate to plot" << endl;
	cout << "    -p (prefix) prefix to candidate files" << endl;
}

int handle_options(int argc,char *argv[])
{
	int c;
    const char *okopts="b:i:l:c:f:t:p:";
    while((c=getopt(argc,argv,okopts))!=-1)
	switch(c)
	{
		case 'i':
			input_eggname=string(optarg);
			break;
		case 'b':
			background_filename=string(optarg);
			break;
		case 'l':
			candidate_list_filename=string(optarg);
			break;
		case 'c':
			convolution_filename=string(optarg);
			break;
		case 'f':
			fbins=atoi(optarg);
			break;
		case 't':
			tbins=atoi(optarg);
			break;
		case 'p':
			prefix=string(optarg);
			break;
	}
	return optind;
}

#include "Waterfall.hh"
#include "Histogram.hh"
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <iostream>
#include <string>
#include <math.h>
#include "correline_utils.hh"
using namespace std;

void print_usage(); //print out instructions
int handle_options(int argc,char *argv[]); //command line options

//----changable parameters----
//int fft_size=1024;
//int fft_size=256;
string input_eggname;
string background_fname;
string convolution_fname;
char format='j';
double frequency_cut_low=20; //in mhz
double frequency_cut_high=80; //in mhz
string prefix="temp_";
double cand_cut=-1;
//----------------------------

Histogram *powers_1;
Histogram *powers_2;
Histogram *powers_combined;
Histogram *convolved_1;
Histogram *convolved_2;
Histogram *convolved_powers;
vector<int> histo_starts;
vector<int> histo_stops;

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
	if(background_fname=="") {
		cerr << "no background file name given" << endl;
		return -1;
	}
	if(convolution_fname=="") {
		cerr << "no convolution file name given" << endl;
		return -1;
	}
		//----load the background----
	Background background;
	background.load(background_fname);
	int fft_size=2*(background.bg_size-1);
	//---------------------------
	//----load the convolution----
	ConvolutionMap cmap;
	cmap.load(convolution_fname);
	//---------------------------
	//----load the egg,prep the correlator-----------
    const Monarch *egg=Monarch::OpenForReading(std::string(input_eggname));
	Correlator correlator;
	correlator.init(egg,fft_size);
	//---------------------------
	//--- make the historgam array
	int n_freq_divisions=5;
	double freq_div_spacing=10;
	powers_1=new Histogram[n_freq_divisions];
	powers_2=new Histogram[n_freq_divisions];
	powers_combined=new Histogram[n_freq_divisions];
	convolved_powers=new Histogram[n_freq_divisions];
	convolved_1=new Histogram[n_freq_divisions];
	convolved_2=new Histogram[n_freq_divisions];
	for(int i=0;i<n_freq_divisions;i++) {
		double divstart=frequency_cut_low+i*freq_div_spacing;
		histo_starts.push_back((divstart*1e6)/correlator.output_waterfall.freq_step);
		double divstop=frequency_cut_low+(i+2)*freq_div_spacing;
		histo_stops.push_back((divstop*1e6)/correlator.output_waterfall.freq_step);
		powers_1[i].init(1000,0,25);
		powers_2[i].init(1000,0,25);
		powers_combined[i].init(1000,0,20);
		convolved_powers[i].init(1000,0,20);
		convolved_1[i].init(1000,0,20);
		convolved_2[i].init(1000,0,20);
	}
	for(int i=0;i<n_freq_divisions;i++) {
		cout << "histo starts " << i << " " << histo_starts[i] << endl;
		cout << "histo stops " << i << " " << histo_stops[i] << endl;
	}
	//------

	//-----figure out where my cuts are----------
	int f_start=int((frequency_cut_low*1e6)/correlator.output_waterfall.freq_step);
	int f_stop=int((frequency_cut_high*1e6)/correlator.output_waterfall.freq_step);
	//---------------------
	

//	int on_event=0;
	double total_time=0;
	double scale=correlator.getScaleFactor();
	int out_size=correlator.output_waterfall.npoints_f;
	
	//----------allocate space--------
	float *channel1_power=new float[out_size*correlator.nffts_per_event];
	float *channel2_power=new float[out_size*correlator.nffts_per_event];
	float *convolution=new float[out_size*correlator.nffts_per_event];
	float *convolution_1=new float[out_size*correlator.nffts_per_event];
	float *convolution_2=new float[out_size*correlator.nffts_per_event];
	PowerCandidateList candidates;
	candidates.max_length=100;
	//
	
	//---------------------
	long total_fft_count=0;
    while(egg->ReadRecord()) {
		//perform power spectrum
		correlator.process_channel1(egg);
		correlator.process_channel2(egg);
		//get power, normalize by background
		for(int i=0;i<correlator.output_waterfall.npoints_t;i++) {
			for(int j=0;j<correlator.output_waterfall.npoints_f;j++) {
				int index=correlator.output_waterfall.getIndex(j,i);
				channel1_power[index]=(correlator.getPowerChannel1(index))*background.background_stdev_invert[j][0];
				channel2_power[index]=(correlator.getPowerChannel2(index))*background.background_stdev_invert[j][1];
				for(int k=0;k<n_freq_divisions;k++) {
					if( (j>histo_starts[k]) && (j<histo_stops[k])) {
						powers_1[k].increment(channel1_power[index]);
						powers_2[k].increment(channel2_power[index]);
						powers_combined[k].increment(0.5*(channel1_power[index]+channel2_power[index]));
					}
				}
			}
		}
		//convolve
		float normalization=1/sqrt((float)cmap.conv_length);
		int tl=cmap.get_t_extent();
		int fl=cmap.get_f_extent();
		for(int i=0;i<correlator.output_waterfall.npoints_t-tl;i++) {
			for(int j=0;j<correlator.output_waterfall.npoints_f-fl;j++) {
				int index=correlator.output_waterfall.getIndex(j,i);
				convolution[index]=0;
				convolution_1[index]=0;
				convolution_2[index]=0;
				for(int k=0;k<cmap.conv_length;k++) {
					int nindex=correlator.output_waterfall.getIndex(j+cmap.conv_fs[k],i+cmap.conv_ts[k]);
					convolution_1[index]+=channel1_power[nindex]+channel2_power[nindex];
					convolution_2[index]+=channel2_power[nindex]+channel2_power[nindex];
					convolution[index]+=0.5*(channel1_power[nindex]+channel2_power[nindex]);
				}
				convolution[index]*=normalization;
				convolution_1[index]*=normalization;
				convolution_1[index]*=normalization;
				//record candidates
				if((j>f_start)&&(j<f_stop))
				if(cand_cut>0)
				if(convolution[index]>cand_cut) {
					double f=correlator.output_waterfall.freq_step*((double)j)/1e6;
					double t=total_time+correlator.output_waterfall.time_step*((double)i);
//					cout << "f t " << f << " " << t << " " << j << " " << i << endl;
					candidates.push(Candidate(f,t,convolution[index]));
				}

				//histogram power
				for(int k=0;k<n_freq_divisions;k++) {
					if( (j>histo_starts[k]) && (j<histo_stops[k])) {
						convolved_powers[k].increment(convolution[index]);
						convolved_1[k].increment(convolution_1[index]);
						convolved_2[k].increment(convolution_2[index]);
					}
				}
			}
		}
//		cout << "timestep: " << ((double)correlator.record_size)t)/(correlator.sampling_rate_mhz*1e6) << endl;

		total_time+=((double)correlator.record_size)/(correlator.sampling_rate_mhz*1e6);
	}
	for(int k=0;k<n_freq_divisions;k++) {
		double divstart=round(frequency_cut_low+k*freq_div_spacing);
		char offnum[256];
		sprintf(offnum,"%d",(int)(divstart));
		string newprefix=prefix+string("_offset")+string(offnum);
		powers_1[k].saveToFile(newprefix+"_power_1_histo.txt");
		powers_2[k].saveToFile(newprefix+"_power_2_histo.txt");
		powers_combined[k].saveToFile(newprefix+"_power_combined_histo.txt");
		convolved_powers[k].saveToFile(newprefix+"_convolved_powers_histo.txt");
		convolved_1[k].saveToFile(newprefix+"_convolved_1_histo.txt");
		convolved_2[k].saveToFile(newprefix+"_convolved_2_histo.txt");
	}
	candidates.saveToFile(prefix+"_candidates.txt");
}

int handle_options(int argc,char *argv[])
{
	int c;
    const char *okopts="i:ab:p:c:d:";
    while((c=getopt(argc,argv,okopts))!=-1)
	switch(c)
	{
		case 'i':
			input_eggname=string(optarg);
			break;
		case 'a':
			format='a';
			break;
		case 'b':
			background_fname=string(optarg);
			break;
		case 'c':
			convolution_fname=string(optarg);
			break;
		case 'p':
			prefix=string(optarg);
			break;
		case 'd':
			cand_cut=atof(optarg);
			break;
		case '?':
			if(index(okopts,optopt)==NULL)
				fprintf(stderr,"{ error: \"unknown option: %c\n, aborting\"}",optopt);
			else
				fprintf(stderr,"{ error: \"option %c does not take an argument, aborting\"}",optopt);
			return -1;

	}
	return optind;
}

void print_usage()
{
	cout << "powerline_electronjager" << endl;
	cout << "" << endl;
	cout << "options: " << endl;
	cout << "	-i (filename) sets the input egg file" << endl;
	cout << "	-a sets format to ascii (default json)" << endl;
	cout << "	-b (filename) sets the background" << endl;
	cout << "   -p (file prefix) sets the prefix to output files" << endl;
}



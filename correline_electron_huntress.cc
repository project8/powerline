#include "Waterfall.hh"
#include "Histogram.hh"
#include <unistd.h>
#include <stdlib.h>
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
double frequency_cut_low=10; //in mhz
double frequency_cut_high=90; //in mhz
string prefix="temp_";
//----------------------------

//------output histograms -------
Histogram correlated_power(400,0,20);
Histogram convolved_power(400,0,20);

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
	//-----figure out where my cuts are----------
	int f_start=int((frequency_cut_low*1e6)/correlator.output_waterfall.freq_step);
	int f_stop=int((frequency_cut_high*1e6)/correlator.output_waterfall.freq_step);
	//---------------------

	//---allocate space----
	Waterfall convolved(correlator.output_waterfall.npoints_t,correlator.output_waterfall.npoints_f,correlator.output_waterfall.freq_step,correlator.output_waterfall.time_step);
	//Waterfall convolved(correlator.output_waterfall.npoints_t,correlator.output_waterfall.npoints_f,correlator.output_waterfall.freq_step,correlator.output_waterfall.time_step);
	//------------------------
	
//	int on_event=0;
//	double total_time=0;
	double scale=correlator.getScaleFactor();
	int out_size=correlator.output_waterfall.npoints_f;
	//---------------------
	long total_fft_count=0;
    while(egg->ReadRecord()) {
		//perform correlation
		correlator.process_event(egg);
		//normalize by the background, fill out power
		for(int i=0;i<correlator.output_waterfall.npoints_t;i++) {
			for(int j=0;j<correlator.output_waterfall.npoints_f;j++) {
	//		for(int j=f_start;j<f_stop;j++) {
				int index=correlator.output_waterfall.getIndex(j,i);
				correlator.output_waterfall.data[index][0]-=background.background_avg[j][0];
				correlator.output_waterfall.data[index][1]-=background.background_avg[j][1];
				correlator.output_waterfall.data[index][0]*=background.background_stdev_invert[j][0];
				correlator.output_waterfall.data[index][1]*=background.background_stdev_invert[j][1];
				//histogram data
				if( (j>=f_start) && (j<=f_stop) )
					correlated_power.increment(correlator.output_waterfall.get_power(index));
			}
		}
		//convolve
		cmap.convolve(&correlator.output_waterfall,&convolved);
		//histogram power
		for(int i=0;i<correlator.output_waterfall.npoints_t;i++) {
			for(int j=f_start;j<f_stop;j++) {
				int index=convolved.getIndex(j,i);
				double p=sqrt(convolved.data[index][0]*convolved.data[index][0]+convolved.data[index][1]*convolved.data[index][1]);
				if( (j>=f_start) && (j<=f_stop) )
					convolved_power.increment(p);
			}
		}
	}
	correlated_power.saveToFile(prefix+"_correlated_power_histogram.txt");
	convolved_power.saveToFile(prefix+"_convolved_power_histogram.txt");
}

int handle_options(int argc,char *argv[])
{
	int c;
    const char *okopts="i:ab:p:c:";
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
	cout << "correline" << endl;
	cout << "returns correlated spectrum" << endl;
	cout << "options: " << endl;
	cout << "	-i (filename) sets the input egg file" << endl;
	cout << "	-a sets format to ascii (default json)" << endl;
	cout << "	-b (filename) sets the background" << endl;
	cout << "   -p (file prefix) sets the prefix to output files" << endl;
}


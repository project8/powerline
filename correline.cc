#include "Waterfall.hh"
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <math.h>
using namespace std;


void print_usage(); //print out instructions
int handle_options(int argc,char *argv[]); //command line options

//----changable parameters----
int fft_size=8192;
//int fft_size=1024;
//int fft_size=256;
//double power_cut_dbm=-40;
//double snr_cut=20;
string input_eggname;
double freq_offset=0;
char format='j';
//----------------------------

//-------storage-----------
double *spectrum_sum_r;
double *spectrum_sum_rr;
double *spectrum_sum_i;
double *spectrum_sum_ii;
//-------------------------

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
    const Monarch *egg=Monarch::OpenForReading(std::string(input_eggname));
	Correlator correlator;
	correlator.init(egg,fft_size);
//	int on_event=0;
//	double total_time=0;
	double scale=correlator.getScaleFactor();

	int out_size=correlator.output_waterfall.npoints_f;
	spectrum_sum_r=new double[out_size];
	spectrum_sum_rr=new double[out_size];
	spectrum_sum_i=new double[out_size];
	spectrum_sum_ii=new double[out_size];
	for(int i=0;i<out_size;i++) {
		spectrum_sum_r[i]=0;
		spectrum_sum_rr[i]=0;
		spectrum_sum_i[i]=0;
		spectrum_sum_ii[i]=0;
	}
	long total_fft_count=0;
    while(egg->ReadRecord()) {
		correlator.process_event(egg);
		for(int i=0;i<correlator.output_waterfall.npoints_t;i++) {
			for(int j=0;j<correlator.output_waterfall.npoints_f;j++) {
				int index=correlator.output_waterfall.getIndex(j,i);
				spectrum_sum_r[j]+=correlator.output_waterfall.data[index][0];
				spectrum_sum_rr[j]+=(correlator.output_waterfall.data[index][0]*correlator.output_waterfall.data[index][0]);
				spectrum_sum_i[j]+=correlator.output_waterfall.data[index][1];
				spectrum_sum_ii[j]+=(correlator.output_waterfall.data[index][1]*correlator.output_waterfall.data[index][1]);
			}
			total_fft_count++;
		}
	}
	if(format=='j') { //json output
		printf("{ \"sampling_rate\": %d , ",correlator.sampling_rate_mhz);
		printf("\"data\": [");
		for(int i=0;i<out_size;i++) {
			if(i!=0) printf(",");
			//this is the power info
			printf("%g",scale*sqrt(spectrum_sum_r[i]*spectrum_sum_r[i]+spectrum_sum_i[i]*spectrum_sum_i[i])/((double)total_fft_count));
			//this is the phase angle info
			//printf("%f",atan2(spectrum_sum_r[i],spectrum_sum_i[i]));
		}
		printf("] }");
	} else { //assume ascii otherwise
		printf("#scale: %g\n",scale);
		printf("#avg_r avg_i stdev_r stdev_i\n");
		for(int i=0;i<out_size;i++) {
			printf("%f ",correlator.output_waterfall.freq_step*((double)i)/1e6);
			double n=((double)total_fft_count);
			double avgr=spectrum_sum_r[i]/n;
			double avgi=spectrum_sum_i[i]/n;
			printf("%g ",scale*avgr);
			printf("%g ",scale*avgi);
			double stdevr=sqrt(spectrum_sum_rr[i]/n-avgr*avgr);
			double stdevi=sqrt(spectrum_sum_ii[i]/n-avgi*avgi);
			printf("%g ",scale*stdevr);
			printf("%g ",scale*stdevi);
			printf("\n");
		}
	}

}

int handle_options(int argc,char *argv[])
{
	int c;
    const char *okopts="i:a";
    while((c=getopt(argc,argv,okopts))!=-1)
	switch(c)
	{
		case 'i':
			input_eggname=string(optarg);
			break;
		case 'a':
			format='a';
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
}

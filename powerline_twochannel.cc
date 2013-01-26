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
//
//
//-------storage-----------
double *channel1_sum;
double *channel2_sum;
double *channel1_sumsq;
double *channel2_sumsq;
double *channel_covariance;
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
	channel1_sum=new double[out_size];
	channel2_sum=new double[out_size];
	channel1_sumsq=new double[out_size];
	channel2_sumsq=new double[out_size];
	channel_covariance=new double[out_size];
	for(int i=0;i<out_size;i++) {
		channel1_sum[i]=0;
		channel2_sum[i]=0;
		channel1_sumsq[i]=0;
		channel2_sumsq[i]=0;
		channel_covariance[i]=0;
	}
	long total_fft_count=0;
    while(egg->ReadRecord()) {
		correlator.process_channel1(egg);
		correlator.process_channel2(egg);
		for(int i=0;i<correlator.output_waterfall.npoints_t;i++) {
			for(int j=0;j<correlator.output_waterfall.npoints_f;j++) {
				int index=correlator.output_waterfall.getIndex(j,i);
				double x1=correlator.fft_output_1[index][0]*correlator.fft_output_1[index][0]+correlator.fft_output_1[index][1]*correlator.fft_output_1[index][1];
				double x2=correlator.fft_output_2[index][0]*correlator.fft_output_2[index][0]+correlator.fft_output_2[index][1]*correlator.fft_output_2[index][1];
				channel1_sum[j]+=x1;
				channel2_sum[j]+=x2;
				channel1_sumsq[j]+=x1*x1;
				channel2_sumsq[j]+=x2*x2;
				channel_covariance[i]+=x1*x2;
			}
			total_fft_count++;
		}
	}
	if(format=='j') { //json output
		printf("{ JSON: \"Not Implemented\"}");
	} else { //assume ascii
		double n=((double)total_fft_count);
		printf("#scale: %g\n",scale);
		printf("#avg_chan1 avg_chan2 stdev_chan1 stdev_chan2\n");
		for(int i=0;i<out_size;i++) {
			printf("%f ",correlator.output_waterfall.freq_step*((double)i)/1e6);
			double avg1=channel1_sum[i]/n;
			double avg2=channel2_sum[i]/n;
			printf("%g ",scale*avg1);
			printf("%g ",scale*avg2);
			double stdev1=sqrt(channel1_sumsq[i]/n-avg1*avg1);
			double stdev2=sqrt(channel2_sumsq[i]/n-avg2*avg2);
			printf("%g %g",scale*stdev1,scale*stdev2);
			double cov=sqrt(channel_covariance[i]/n-avg1*avg2);
			printf(" %g\n",scale*cov);
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
	cout << "powerline_twochannel" << endl;
	cout << "returns power spectrum and standard deviations for both channels" << endl;
	cout << "options: " << endl;
	cout << "	-i (filename) sets the input egg file" << endl;
	cout << "	-a sets format to ascii (default json)" << endl;
}

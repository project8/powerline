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
//int fft_size=1024;
int fft_size=256;
double power_cut_dbm=-16;
string input_eggname;
double freq_offset=0;
char format='j';
//----------------------------

//----structures I fill out---
double *uncorrelated_shape_1=NULL;
double *uncorrelated_shape_2=NULL;
double *correlated_average=NULL;
double *correlated_average_counter=NULL;
double *avg_1=NULL;
double *avg_2=NULL;
//----------------------------

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
	int on_event=0;
	double total_time=0;
	double scale=correlator.getScaleFactor();
	double correlation_cut=pow(pow10(power_cut_dbm/10.0)/scale,2.0);

	int out_size=fft_size/2+1;
	uncorrelated_shape_1=new double[out_size];
	uncorrelated_shape_2=new double[out_size];
	correlated_average=new double[out_size];
	correlated_average_counter=new double[out_size];
	avg_1=new double[out_size];
	avg_2=new double[out_size];
	long n_uncorrelated_shape_avgs=0;
	for(int i=0;i<out_size;i++) {
		uncorrelated_shape_1[i]=0;
		uncorrelated_shape_2[i]=0;
		correlated_average[i]=0;
		correlated_average_counter[i]=0;
		avg_1[i]=0;
		avg_2[i]=0;
	}
	int sampling_rate_mhz=correlator.sampling_rate_mhz;
    while(egg->ReadRecord()) {
		correlator.process_event(egg);
		for(int i=0;i<correlator.output_waterfall.npoints_t;i++) {
			double max_power=0;
			int max_power_pos=0;
			for(int j=0;j<correlator.output_waterfall.npoints_f;j++) {
				int index=correlator.output_waterfall.getIndex(j,i);
//				fftwf_complex a=correlator.fft_output_1[index]
//				fftwf_complex b=correlator.fft_output_2[index]
				double power=correlator.output_waterfall.get_power_squared(index);
				if(power>max_power) {
					max_power=power;
					max_power_pos=j;
				}
				/*
				if(power>correlation_cut)
					cout << total_time+correlator.output_waterfall.time_step*((double)i) << " " << correlator.output_waterfall.freq_step*((double)j) << " " << sqrt(power)*scale << endl;
					*/

			}
			//if the sweeper is turned off, use this to measure uncorrelated power
			if(max_power<correlation_cut) {
			for(int j=0;j<correlator.output_waterfall.npoints_f;j++) 
			{
				int index=correlator.output_waterfall.getIndex(j,i);
				double p_a=correlator.getPowerChannel1(index);
				double p_b=correlator.getPowerChannel2(index);
				uncorrelated_shape_1[j]+=p_a;
				uncorrelated_shape_2[j]+=p_b;
			}
			n_uncorrelated_shape_avgs++;
			} else {
			//otherwise, average the correlated power
			//a crude handling of when the sweep is between bins here
			if(max_power_pos==0) max_power_pos=1;
			if(max_power_pos==out_size-1) max_power_pos=out_size-2;
			int aindex=correlator.output_waterfall.getIndex(max_power_pos-1,i);
			int bindex=correlator.output_waterfall.getIndex(max_power_pos,i);
			int cindex=correlator.output_waterfall.getIndex(max_power_pos+1,i);
			double a=sqrt(correlator.output_waterfall.get_power_squared(aindex));
			double b=sqrt(max_power);
			double c=sqrt(correlator.output_waterfall.get_power_squared(cindex));
			double sum=a+b+c;
			double afrac=a/sum;
			double bfrac=b/sum;
			double cfrac=c/sum;
			avg_1[max_power_pos-1]+=afrac*correlator.getPowerChannel1(aindex);
			avg_2[max_power_pos-1]+=afrac*correlator.getPowerChannel2(aindex);
			correlated_average[max_power_pos-1]+=a*afrac;
			correlated_average_counter[max_power_pos-1]+=afrac;
			avg_1[max_power_pos]+=bfrac*correlator.getPowerChannel1(bindex);
			avg_2[max_power_pos]+=bfrac*correlator.getPowerChannel2(bindex);


			correlated_average[max_power_pos]+=b*bfrac;
			correlated_average_counter[max_power_pos]+=bfrac;
			avg_1[max_power_pos+1]+=cfrac*correlator.getPowerChannel1(cindex);
			avg_2[max_power_pos+1]+=cfrac*correlator.getPowerChannel2(cindex);

			correlated_average[max_power_pos+1]+=c*cfrac;
			correlated_average_counter[max_power_pos+1]+=cfrac;
			}
		}
		total_time+=correlator.getEventDuration();
	}
	//print out json
//	/*
	if(format=='j') {
		printf("{ \"sampling_rate\": %d , ",sampling_rate_mhz);
		printf("\"data\": [");
		for(int i=0;i<out_size;i++) {
			if(i!=0) printf(",");
			if(correlated_average_counter[i]!=0)
			//printf("%f",scale*correlated_average[i]/correlated_average_counter[i]);
			printf("%f",scale*avg_1[i]/correlated_average_counter[i]);
			else
			printf("0");
		}
		printf("] }");
	} else { //assume ascii otherwise
		for(int i=0;i<(fft_size/2+1);i++) {
			cout << freq_offset+(correlator.output_waterfall.freq_step*((double)i))/1e6 << " ";
			if(correlated_average_counter[i]>0) {
				cout << scale*(correlated_average[i]/(correlated_average_counter[i])) << " ";
				cout << scale*avg_1[i]/(correlated_average_counter[i]) << " ";
				cout << scale*avg_2[i]/(correlated_average_counter[i]) << " ";
			} else
				cout << "0 0 0 ";
			cout << scale*uncorrelated_shape_1[i]/((double)n_uncorrelated_shape_avgs) << " ";
			cout << scale*uncorrelated_shape_2[i]/((double)n_uncorrelated_shape_avgs) << endl;
		}

	}
//	*/
	
}


	
int handle_options(int argc,char *argv[])
{
	int c;
    const char *okopts="i:o:a";
    while((c=getopt(argc,argv,okopts))!=-1)
	switch(c)
	{
		case 'i':
			input_eggname=string(optarg);
			break;
		case 'o':
			freq_offset=atof(optarg);
			break;
		case 'a':
			format='a';
			break;
		case '?':
			//FILE *errfl=fopen("errfile.log","w");
			if(index(okopts,optopt)==NULL)
				fprintf(stderr,"{ error: \"unknown option: %c\n, aborting\"}",optopt);
			else
				fprintf(stderr,"{ error: \"option %c does not take an argument, aborting\"}",optopt);
			//fclose(errfl);
			return -1;

	}
	return optind;
}

void print_usage()
{
	cout << "dpph_search" << endl;
	cout << "experimental dpph search" << endl;
	cout << "options: " << endl;
	cout << "	-i (filename) sets the input egg file" << endl;
	cout << "	-a sets format to ascii (default json)" << endl;
}

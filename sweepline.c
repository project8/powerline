#include <fftw3.h>
#include <unistd.h>
#include "Monarch.hpp"
//#include "monarch.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


/*---Configurable Settings---*/
int fft_size=1024;
int max_number_of_events=1024;
//double sweep_cut_dbm=-40;
double sweep_cut_dbm=-15;
char format='j'; //b=binary, a=ascii, j=json
char eggname[512];
int on_channel=1;
//int on_channel=2;
/*---------------------------*/

/*---FFT buffers and such--*/
int nffts_per_event;
float *fft_input;
fftwf_plan fft_plan=NULL;
fftwf_complex *fft_output;
int fft_output_size;
double *output_powerspectrum;
int *output_powercount;
int sampling_rate_mhz;
/*----------------*/


void print_usage(); //print out instructions
int handle_options(int argc,char *argv[]); //command line options

int main(int argc,char *argv[]) 
{
	//this is C, so I have to declare all sorts of variables in advance
	int i;
	int j;
	eggname[0]='\0';

	//handle the command line options
	int onindex;
    if(((onindex=handle_options(argc,argv))==-1)||(argc-onindex<0))
		{print_usage(); return -1;};
//	char *eggname=argv[onindex];
	if(eggname[0]=='\0') {
		fprintf(stderr,"no input file given, use -i option\n");
		return -1;
	}
	double sweep_cut_mw=pow(10.0,sweep_cut_dbm*0.1);

	//open the egg
/*
	struct egg current;
	mBreakEgg(eggname,&current);
	mParseEggHeader(&current);
	sampling_rate_mhz=current.data->sample_rate;
*/
	//Monarch *egg=Monarch::Open(std::string(eggname),ReadMode);
	const Monarch *egg=Monarch::OpenForReading(std::string(eggname));
	egg->ReadHeader();

	const MonarchHeader *eggheader=egg->GetHeader();
	const MonarchRecord *event;
	sampling_rate_mhz=eggheader->GetAcqRate();


	//decide the optimal size for ffts and allocate memory
//	printf("record size: %d\n",eggheader->GetRecordSize());
	if(eggheader->GetRecordSize()<(unsigned int)fft_size) {
	//if(current.data->record_size<fft_size) {
		fprintf(stderr,"fft size larger than record size.  make fft size smaller. aborting");
		return -1;
	}
	//nffts_per_event=current.data->record_size/fft_size;
	nffts_per_event=eggheader->GetRecordSize()/fft_size;
	fft_output_size=fft_size/2+1;
	//fft_input=fftwf_alloc_real(fft_size*nffts_per_event);
	fft_input=(float*)fftwf_malloc(sizeof(float)*nffts_per_event*fft_size);
	//fft_output=fftwf_alloc_complex(fft_size*nffts_per_event);
	fft_output=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	output_powerspectrum=(double*)malloc(sizeof(double)*fft_output_size);
	output_powercount=(int*)malloc(sizeof(int)*fft_output_size);
	for(i=0;i<fft_output_size;i++) {
		output_powerspectrum[i]=0;
		output_powercount[i]=0;
	}

	//normalize to power in milliwatts
	//1000 (milliwatts/watt) * 0.5 (volts fullscale)/256 (levels) / (sqrt(fft_length)*(number of averages) / 50 ohms
	double normalization=2.0*(1000.0*0.5*0.5/(256.0*256.0))*(1.0/(((double)fft_size*fft_size)))/50.0;

	//create the fft plan
	fft_plan=fftwf_plan_many_dft_r2c(1,&fft_size,nffts_per_event,fft_input,NULL,1,fft_size,fft_output,NULL,1,fft_output_size,FFTW_ESTIMATE);

	//perform ffts
	int on_event=0;
	int nffts_so_far=0;
	double power=0;
	//while((mHatchNextEvent(&current)!=1)&&(on_event<=max_number_of_events)) {

//	while((event=egg->GetNextEvent())!=NULL&&(on_event<=max_number_of_events)) {
	while(egg->ReadRecord()) {
		if(on_channel==1) {
			event=egg->GetRecordOne();
		} else {
			event=egg->GetRecordTwo();
		}

		//convert data to floats
		//for(i=0;i<current.data->record_size;i++)
		for(i=0;i<eggheader->GetRecordSize();i++)
		{
			//fft_input[i]=(float)(current.data->record[i])-128.0;
			fft_input[i]=(float)(event->fDataPtr[i])-128.0;
		}
		//perform the ffts
		fftwf_execute(fft_plan);
		//pack in to power spectrum
		int on_pt=0;
		for(i=0;i<nffts_per_event;i++)
		for(j=0;j<fft_output_size;j++) {
			power=normalization*(fft_output[on_pt][0]*fft_output[on_pt][0]+fft_output[on_pt][1]*fft_output[on_pt][1]);
			if(power>sweep_cut_mw) {
				output_powerspectrum[j]+=power;
				output_powercount[j]++;
			}
			on_pt++;
		}
		nffts_so_far+=nffts_per_event;
	}
	//normalize each bin by the number of spectra that contibuted to it
	for(i=0;i<fft_output_size;i++)  {
		if(output_powercount[i]!=0)
		output_powerspectrum[i]/=((double)output_powercount[i]);
	}

	//print out result
	if(format=='j') { //ASCII output, JSON
		printf("{ \"sampling_rate\": %d , ",sampling_rate_mhz);
		printf("\"data\": [");
		for(i=0;i<fft_output_size;i++) {
			if(i!=0) printf(",");
			printf("%g",output_powerspectrum[i]);
		} 
		printf("] }");
	} else if(format=='a') { 
		for(i=0;i<fft_output_size;i++)
			printf("%g %g\n",sampling_rate_mhz*((double)i)/((double)(2*fft_output_size)),output_powerspectrum[i]);
	} else { //binary
		fwrite(output_powerspectrum,sizeof(double),fft_output_size,stdout);
	}

	//clean up
	fftwf_destroy_plan(fft_plan);
	fftwf_free(fft_input);
	fftwf_free(fft_output);
	free(output_powerspectrum);
	free(output_powercount);
	//mCleanUp(&current);
	return 0;
}

void print_usage()
{
	printf("sweepline\n");
	printf("prints out a power spectrum from an egg file");
	printf("but only includes power over a certain value, for use with sweepers\n");
	printf("Usage: sweepline [options]\n");
	printf("  options:\n");
	printf("  -i (filename) sets the input egg file  MANDATORY\n");
	printf("  -a     sets output to plain ASCII (default JSON)\n");
	printf("  -b     sets output to binary (default JSON)\n");
	printf("  -f (integer)  sets the number of points in the fft (default %d)\n",fft_size);
	printf("  -n (integer)  sets the maximum number of events to scan (default %d)\n",max_number_of_events);
	printf("  -c (1 or 2)  sets the channed being used (default %d)\n",on_channel);
	printf("  -g (floating point) sets the cut in dBm below which power will not be averaged (default %g)\n",sweep_cut_dbm);
}

int handle_options(int argc,char *argv[])
{
    int c;
    const char *okopts="abf:n:i:c:g:";
    while((c=getopt(argc,argv,okopts))!=-1)
	switch(c)
	{
		case 'a':
			format='a';
			break;
		case 'b':
			format='b';
			break;
		case 'g':
			sweep_cut_dbm=atof(optarg);
			break;
		case 'c':
			on_channel=atoi(optarg);
			break;
		case 'f':
			fft_size=atoi(optarg);
			if((fft_size<2)||(fft_size>1024*1024)) {
				fprintf(stderr,"fft size is insane.  aborting.\n");
				return -1;
			}
			break;
		case 'n':
			max_number_of_events=atoi(optarg);
			if(max_number_of_events<1) {
				fprintf(stderr,"max number of events is insane. aborting.\n");
				return -1;
			}
			break;
		case 'i':
			strcpy(eggname,optarg);
			break;
		case '?':
			if(index(okopts,optopt)==NULL)
				fprintf(stderr,"unknown option: %c\n, aborting",optopt);
			else
				fprintf(stderr,"option %c does not take an argument, aborting\n",optopt);
			return -1;
	}
	return optind;
}

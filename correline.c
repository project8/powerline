#include <fftw3.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Monarch.hpp"

/*---Configurable Settings---*/
int fft_size=1024;
int max_number_of_events=1024;
char format='j'; //b=binary, a=ascii, j=json
char eggname[512];
/*---------------------------*/

/*---FFT buffers and such--*/
int nffts_per_event;
float *fft_input;
fftwf_plan fft_plan_1=NULL;
fftwf_plan fft_plan_2=NULL;
fftwf_complex *fft_output_1;
fftwf_complex *fft_output_2;
fftwf_complex *correlation_average;
fftwf_complex *correlation_squared;
int fft_output_size;
int sampling_rate_mhz;
/*----------------*/


void print_usage(); //print out instructions
int handle_options(int argc,char *argv[]); //command line options

int main(int argc,char *argv[]) 
{
	//this is C, so I have to declare all sorts of variables in advance
	//update, I guess it's C++ now so this is a luxury
	int i;
	int j;
	double stdev1;
	double stdev2;
	eggname[0]='\0';

	//handle the command line options
	int onindex;
    if(((onindex=handle_options(argc,argv))==-1)||(argc-onindex<0))
		{print_usage(); return -1;};
	if(eggname[0]=='\0') {
		fprintf(stderr,"no input file given, use -i option\n");
		return -1;
	}

	//open the egg
	const Monarch *egg=Monarch::OpenForReading(std::string(eggname));
	egg->ReadHeader();
	const MonarchHeader *eggheader=egg->GetHeader();
	const MonarchRecord *event;
	sampling_rate_mhz=eggheader->GetAcqRate();

	int record_size=eggheader->GetRecordSize()/2;

	//decide the optimal size for ffts and allocate memory
	if(record_size<fft_size) {
		fprintf(stderr,"fft size larger than record size.  make fft size smaller. aborting");
		return -1;
	}
	nffts_per_event=record_size/fft_size;
	fft_output_size=fft_size/2+1;
	fft_input=(float*)fftwf_malloc(sizeof(float)*nffts_per_event*fft_size);
	fft_output_1=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	fft_output_2=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	correlation_average=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	correlation_squared=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	for(i=0;i<fft_output_size;i++)	{ 
		correlation_average[i][0]=0;
		correlation_average[i][1]=0;
		correlation_squared[i][1]=0;
		correlation_squared[i][1]=0;
		}
	//output_powerspectrum=(double*)malloc(sizeof(double)*fft_output_size);
	//for(i=0;i<fft_output_size;i++) output_powerspectrum[i]=0;

	//create the fft plan
	fft_plan_1=fftwf_plan_many_dft_r2c(1,&fft_size,nffts_per_event,fft_input,NULL,1,fft_size,fft_output_1,NULL,1,fft_output_size,FFTW_ESTIMATE);
	fft_plan_2=fftwf_plan_many_dft_r2c(1,&fft_size,nffts_per_event,fft_input,NULL,1,fft_size,fft_output_2,NULL,1,fft_output_size,FFTW_ESTIMATE);

	//perform ffts
	//int last_channel=-1; //-1 means last channel buffer empty
	int nffts_so_far=0;
	//while((mHatchNextEvent(&current)!=1)&&(on_event<=max_number_of_events)) {
	//while((event=egg->GetNextEvent())!=NULL&&(on_event<=max_number_of_events)) {
	while(egg->ReadRecord()) {
		event=egg->GetRecordOne();
		//convert data to floats
		for(i=0;i<record_size;i++)
			fft_input[i]=(float)(event->fDataPtr[i])-128.0;
		//perform the ffts
		fftwf_execute(fft_plan_1);
		//and the other channel
		event=egg->GetRecordTwo();
		for(i=0;i<record_size;i++)
			fft_input[i]=(float)(event->fDataPtr[i])-128.0;
		fftwf_execute(fft_plan_2);
		int on_pt=0;
		for(i=0;i<nffts_per_event;i++)
		for(j=0;j<fft_output_size;j++) {
			fftwf_complex mult;
			mult[0]=fft_output_1[on_pt][0]*fft_output_2[on_pt][0]+fft_output_1[on_pt][1]*fft_output_2[on_pt][1];
			mult[1]=fft_output_1[on_pt][0]*fft_output_2[on_pt][1]-fft_output_1[on_pt][1]*fft_output_2[on_pt][0];
			correlation_average[j][0]+=mult[0];
			correlation_average[j][1]+=mult[1];
			correlation_squared[j][0]+=mult[0]*mult[0];
			correlation_squared[j][1]+=mult[1]*mult[1];
			on_pt++;
		}
		nffts_so_far+=nffts_per_event;
	}
	//1000 (milliwatts/watt) * 0.5 (volts fullscale)/256 (levels) / (sqrt(fft_length)*(number of averages) / 50 ohms
	double normalization=2.0*(1000.0*0.5*0.5/(256.0*256.0))*(1.0/(((double)fft_size*fft_size)*((double)nffts_so_far)))/50.0;
	double square_norm=normalization*normalization*((double)nffts_so_far);
	for(i=0;i<fft_output_size;i++) {
		correlation_average[i][0]*=normalization;
		correlation_average[i][1]*=normalization;
		correlation_squared[i][0]*=square_norm;
		correlation_squared[i][1]*=square_norm;
	}


	//print out result
	if(format=='j') { //ASCII output, JSON
		printf("{ \"sampling_rate\": %d , ",sampling_rate_mhz);
		printf("\"data\": [");
		for(i=0;i<fft_output_size;i++) {
			if(i!=0) printf(",");
			printf(" [ %g , %g ] ",correlation_average[i][0],correlation_average[i][1]);
		} 
		printf("] }");
	} else if(format=='a') { 
		printf("#chan1_correl chan2_correl chan1_stdev chan2_stdev\n");
		for(i=0;i<fft_output_size;i++) {
			//printf("%g %g %g\n",sampling_rate_mhz*((double)i)/((double)(2*fft_output_size)),correlation_average[i][0],correlation_average[i][1]);
			stdev1=sqrt(correlation_squared[i][0]-correlation_average[i][0]*correlation_average[i][0]);
			stdev2=sqrt(correlation_squared[i][1]-correlation_average[i][1]*correlation_average[i][1]);
			printf("%g %g %g %g %g\n",sampling_rate_mhz*((double)i)/((double)(2*fft_output_size)),correlation_average[i][0],correlation_average[i][1],stdev1,stdev2);
		}
	} else { //binary
		fwrite(correlation_average,sizeof(fftwf_complex),fft_output_size,stdout);
	}

	//clean up
	egg->Close();
	fftwf_destroy_plan(fft_plan_1);
	fftwf_destroy_plan(fft_plan_2);
	fftwf_free(fft_input);
	fftwf_free(fft_output_1);
	fftwf_free(fft_output_2);
	fftwf_free(correlation_average);
	//mCleanUp(&current);
	return 0;
}

void print_usage()
{
	printf("correline\n");
	printf("prints out a correlation spectrum from an egg file");
	printf("Usage: powerline [options]\n");
	printf("  options:\n");
	printf("  -i (filename) sets the input egg file  MANDATORY\n");
	printf("  -a     sets output to plain ASCII (default JSON)\n");
	printf("  -b     sets output to binary (default JSON)\n");
	printf("  -f (integer)  sets the number of points in the fft (default %d)\n",fft_size);
	printf("  -n (integer)  sets the maximum number of events to scan (default %d)\n",max_number_of_events);
}

int handle_options(int argc,char *argv[])
{
    int c;
    const char *okopts="abf:n:i:";
    while((c=getopt(argc,argv,okopts))!=-1)
	switch(c)
	{
		case 'a':
			format='a';
			break;
		case 'b':
			format='b';
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

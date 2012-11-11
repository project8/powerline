#include <fftw3.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
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


	//decide the optimal size for ffts and allocate memory
	if(eggheader->GetRecordSize()<(unsigned int)fft_size) {
		fprintf(stderr,"fft size larger than record size.  make fft size smaller. aborting");
		return -1;
	}
	nffts_per_event=eggheader->GetRecordSize()/fft_size;
	fft_output_size=fft_size/2+1;
	fft_input=(float*)fftwf_malloc(sizeof(float)*nffts_per_event*fft_size);
	fft_output_1=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	fft_output_2=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	correlation_average=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	for(i=0;i<fft_output_size;i++)	{ 
		correlation_average[i][0]=0;
		correlation_average[i][1]=0;
		}
	//output_powerspectrum=(double*)malloc(sizeof(double)*fft_output_size);
	//for(i=0;i<fft_output_size;i++) output_powerspectrum[i]=0;

	//create the fft plan
	fft_plan_1=fftwf_plan_many_dft_r2c(1,&fft_size,nffts_per_event,fft_input,NULL,1,fft_size,fft_output_1,NULL,1,fft_output_size,FFTW_ESTIMATE);
	fft_plan_2=fftwf_plan_many_dft_r2c(1,&fft_size,nffts_per_event,fft_input,NULL,1,fft_size,fft_output_2,NULL,1,fft_output_size,FFTW_ESTIMATE);

	//perform ffts
	int last_channel=-1; //-1 means last channel buffer empty
	int nffts_so_far=0;
	//while((mHatchNextEvent(&current)!=1)&&(on_event<=max_number_of_events)) {
	//while((event=egg->GetNextEvent())!=NULL&&(on_event<=max_number_of_events)) {
	while(egg->ReadRecord()) {
		event=egg->GetRecordOne();
		//convert data to floats
		for(i=0;i<eggheader->GetRecordSize();i++)
			fft_input[i]=(float)(event->fDataPtr[i])-128.0;
		//perform the ffts
		fftwf_execute(fft_plan_1);
		//and the other channel
		event=egg->GetRecordTwo();
		for(i=0;i<eggheader->GetRecordSize();i++)
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
			on_pt++;
		}
		nffts_so_far+=nffts_per_event;
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
		for(i=0;i<fft_output_size;i++)
			printf("%g %g %g\n",sampling_rate_mhz*((double)i)/((double)(2*fft_output_size)),correlation_average[i][0],correlation_average[i][1]);
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

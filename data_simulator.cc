#include <strings.h>
#include <unistd.h>
#include "Monarch.hpp"
#include <string>
#include <math.h>
#include <stdlib.h>
#include "correline_utils.hh"
using namespace std;

double spare_grand;
bool has_spare_grand=false;
bool empty=false;
double getGaussianRand(double mean,double sigma);
double getUniformRand(double range);
int handle_options(int argc,char *argv[]);
string outfname="test.egg";
string truthfname="truth.txt";
PowerCandidateList candidates;
int nrecords=6;
int signalevery=1;
double snr=2.0;
double energy_loss=1e-15; //watts
double wattstoevs=(1/1.6e-19); // eV/s per Watt
double me=511000; //eV
double ee=32000; //eV
double f0=26e9; //Hz
double power=1e-16; //Watts

void print_usage(); //print out instructions

int main(int argc,char *argv[])
{
    double freq_per_second=f0*(energy_loss*wattstoevs)/(me+ee); //Hz/s
	cout << "frequency change per second " << freq_per_second << endl;
	srand(time(NULL));
    //handle the command line options
    int onindex;
    if(((onindex=handle_options(argc,argv))==-1)||(argc-onindex<0))
	{cerr << "argument error" << endl;
		return -1;};

	//int record_size=262144;
	int record_size=4194304;
	int event_length=40000; //200e-6 s * 200 MHz
	//int event_length=1000000; //a long time
	int sampling_rate=200e6;
/*
if each sample has a gaussian with stdev sigma
and N are summed together * sin(omega t)
the expected value is 0 with stdev 0.5*sigma*sqrt(N)
so the average value squared is N*sigma^2/4

if each sample has A*sin(omega*t) 
and N are summed together * sin(omega*t)
then the expected value is 0.5*N*A
so the average value squared is N^2*A^2/4

given that SNR is defined as P_sig/(KT) * sqrt (t/b)
b=1/t
N=r*t where r is sampling rate
P_sig/(KT) *N/r

compare with above
SNR=N*(A^2/sigma^2)

So I choose KT/r = 1.4e-23 (J/K) * (100K) * (2e8 Hz) = 2e-13 Watts
if I expect my signal to have 1e-15 Watts, then
A^2/sigma^2=5e-3
if sigma=1
so A=0.07

 */
	double noise_magnitude=pow(2,6);
	Monarch *egg=Monarch::OpenForWriting(outfname);
	MonarchHeader *header=egg->GetHeader();
	header->SetAcqRate(sampling_rate/1e6);
	header->SetRecordSize(record_size);
	header->SetAcqTime(100);
	header->SetAcqMode(sTwoChannel);
	if(!egg->WriteHeader()) {
		cerr << "failed to write header" << endl;
	}
	double noise_power=100*1.4e-23*(sampling_rate/2);
	snr=sqrt(power/noise_power);
	cout << "snr is " << snr << endl;

	double phase_diff=2*M_PI*getUniformRand(1.0);
	for(int onrecord=0;onrecord<nrecords;onrecord++) {
		int event_start=getUniformRand(record_size-event_length);
		double event_freq=(20+60*getUniformRand(1.0))*1e6;
		candidates.push(Candidate(event_freq/1e6,(((double)onrecord)*((double)record_size)+((double)event_start))/((double)sampling_rate),1.0));
//		cout << "event at " << (((double)onrecord)*((double)record_size)+((double)event_start))/(1e6*(double)sampling_rate) << " s" << endl;
//		const MonarchRecord *r1=egg->GetRecordOne();
//		const MonarchRecord *r2=egg->GetRecordTwo();
		MonarchRecord *r=egg->GetRecordInterleaved();
		for(int i=0;i<record_size;i++) {
			double x1=getGaussianRand(128.0,noise_magnitude);
			double x2=getGaussianRand(128.0,noise_magnitude);
			if(!empty)
			if(onrecord%signalevery==0)
			if((i>=event_start)&&(i<event_start+event_length)) {
			    double x=((double)(i-event_start))/sampling_rate;
			    double phase=x*2*M_PI*(event_freq+x*freq_per_second);
				x1+=snr*noise_magnitude*sin(phase);
				x2+=snr*noise_magnitude*cos(phase+phase_diff);
			}
			//clipping
			if(x1>256) x1=256;
			if(x1<0) x1=0;
			if(x2>256) x2=256;
			if(x2<0) x2=0;
			r->fDataPtr[2*i]=(unsigned char)x1;
			r->fDataPtr[2*i+1]=(unsigned char)x2;
			//r1->fDataPtr[i]=(unsigned char)x1;
			//r2->fDataPtr[i]=(unsigned char)x2;
		}
		if(!egg->WriteRecord()) {
			cerr << "failed to write record" << endl;
		}
	}
	egg->Close();

	candidates.saveToFile(truthfname);


}

double getUniformRand(double range)
{
    return range*((double)rand())/((double)RAND_MAX);
}

double getGaussianRand(double mean,double sigma)
{
    //Polar form of the Box-Muller transformation
    //taken from http://www.taygeta.com/random/gaussian.html
    if(has_spare_grand) {has_spare_grand=false; return spare_grand*sigma+mean;}
   double x1, x2, w, y1;
   do {
             x1 = 2.0 * getUniformRand(1.0) - 1.0;
             x2 = 2.0 * getUniformRand(1.0) - 1.0;
             w = x1 * x1 + x2 * x2;
   } while ( w >= 1.0 );
   w = sqrt( (-2.0 * log( w ) ) / w );
   y1 = x1 * w;
    spare_grand= x2 * w;
    has_spare_grand=true;
   return y1*sigma+mean;
}

void print_usage()
{
	cout << "project 8 data simulator" << endl;
	cout << "makes a fake egg file" << endl;
	cout << "usage data_simulator [options]" << endl;
	cout << "-e data file will be empty of candidates regardless of other options" << endl;
	cout << "-n (number) this many events will be produced" << endl;
	cout << "-o (filename) this is the output egg file name" << endl;
	cout << "-t (filename) this is where the mc truth is stored, locations of candidates" << endl;
	cout << "-p (number) this is the power of the signals in watts" << endl;
	cout << "-m (number) signals are made every m events" << endl;
}

int handle_options(int argc,char *argv[])
{
    int c;
    const char *okopts="eo:t:p:n:m:";
    while((c=getopt(argc,argv,okopts))!=-1)
	switch(c)
	{
		case 'e':
			empty=true;
			break;
		case 'o':
			outfname=string(optarg);
			break;
		case 't':
			truthfname=string(optarg);
			break;
		case 'n':
			nrecords=atoi(optarg);
			break;
		case 'p':
			power=atof(optarg);
			break;
		case 's':
			snr=atof(optarg);
			break;
		case 'm':
			signalevery=atoi(optarg);
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


#include <fftw3.h>
#include "Monarch.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

//Contains information from a series of limited time FFTs
class Waterfall
{
public:
	Waterfall() {data=NULL;};
	Waterfall(int npt,int npf,double fstep,double tstep);
	~Waterfall();
	void init(int npt,int npf,double fstep,double tstep);

	int getIndex(int nfreq,int ntime) {return ntime*npoints_f+nfreq;};
	int getFrequencyBin(double freq) {
		return freq/freq_step;};
	int getTimeBin(double time) {
		return time/time_step;};

	void saveSubWaterfall(int fstart,int fend,int tstart,int tend,double time_offset,string fname);

	void subtractBackground(fftwf_complex *spectrum);
	void multiplyByFrequency(fftwf_complex *spectrum);

	double get_power_squared(int index);


	fftwf_complex *data;

	int npoints_f;     //number of frequency points
	int npoints_t;     //number of time points
	double freq_step;  //frequency bin size
	double time_step;  //time bin size
};

//Not really a class, just a collection of stuff to make producing
//a correlation waterfall plot convenient
class Correlator
{
public:
	void init(const Monarch *egg,int myfftsize);
	void process_event(const Monarch *egg);
	double getScaleFactor(); //returns units of mW

	float *fft_input;
	fftwf_plan fft_plan_1;
	fftwf_plan fft_plan_2;
	fftwf_complex *fft_output_1;
	fftwf_complex *fft_output_2;
	Waterfall output_waterfall;

	double getPowerChannel1(int index);
	double getPowerChannel2(int index);

	int fft_size;
	int nffts_per_event;
	int fft_output_size;
	int record_size;
	int sampling_rate_mhz; //sampling rate in MHz
	
	double getEventDuration() {
		return ((double)(nffts_per_event*fft_size))/(((double)sampling_rate_mhz)*1e6);
	}

};
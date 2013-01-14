#pragma once
#include <fftw3.h>
#include <string>
#include <vector>
#pragma once
#include <set>
#include <iostream>
#include <sstream>
#include <fstream>
#include "Waterfall.hh"
using namespace std;

class Background
{
public:
	~Background();
	bool load(string fname);

	fftwf_complex *background_avg;
	fftwf_complex *background_stdev;
	fftwf_complex *background_stdev_invert;
	int bg_size;
};

class ConvolutionMap
{
public:
	int conv_ts[100];
	int conv_fs[100];
	int conv_length;
	int get_t_extent();
	int get_f_extent();
	void convolve(Waterfall *input,Waterfall *output);
	bool load(string fname);
};

class Candidate
{
public:
	Candidate() {};
	Candidate(double f,double t,double m) {
		frequency=f;
		time=t;
		magnitude=m;
	};
	double frequency;
	double time;
	double magnitude;
	int event_no;
	//--optional parameters--
	double probability; //probability of being a signal

	int getEventNo(int sampling_rate_mhz,int event_size);
};

class SortCandidatesByPower {
public:
	bool operator()(const Candidate &a,const Candidate &b);
};

class PowerCandidateList : public set<Candidate,SortCandidatesByPower>
{
public:
	PowerCandidateList() {max_length=1000000;};
	void push(const Candidate &c);
	unsigned int max_length;
	void saveToFile(string fname,long long npoints=0);
	PowerCandidateList condenseByDistance(double freqscale,double timescale);
};

vector<Candidate> loadCandidateFile(string fname);
vector<Candidate> loadCandidateFile(string fname,long long &npoints);

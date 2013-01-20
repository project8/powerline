#pragma once
#include <string>
using namespace std;

class Histogram
{
public:
	Histogram() {data=NULL; size=0;};
    Histogram(int sz,float mmin,float mmax);
    ~Histogram();
	void init(int sz,float mmin,float mmax);
    void saveToFile(string fname);

    void increment(float x);
    double getx(int bin);
	double getbin(double x);
    double calc_mean();
    double calc_sigma();

    double *data;
    int size;
    float min;
    float max;
};

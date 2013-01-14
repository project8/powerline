#include "Waterfall.hh"
#include <math.h>

Waterfall::Waterfall(int npt,int npf,double fstep,double tstep)
{
	data=NULL;
	init(npt,npf,fstep,tstep);
}
	
void Waterfall::init(int npt,int npf,double fstep,double tstep)
{
	npoints_t=npt;
	npoints_f=npf;
	freq_step=fstep;
	time_step=tstep;
    data=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*npoints_t*npoints_f);
}
	
Waterfall::~Waterfall()
{
	fftwf_free(data);
}

void Waterfall::subtractBackground(fftwf_complex *spectrum)
{
	int on_pt=0;
	for(int i=0;i<npoints_t;i++)
   	for(int j=0;j<npoints_f;j++) {
		data[on_pt][0]-=spectrum[j][0];
		data[on_pt][1]-=spectrum[j][1];
		++on_pt;
	}
}

void Waterfall::multiplyByFrequency(fftwf_complex *spectrum)
{
	int on_pt=0;
	for(int i=0;i<npoints_t;i++)
   	for(int j=0;j<npoints_f;j++) {
		data[on_pt][0]*=spectrum[j][0];
		data[on_pt][1]*=spectrum[j][1];
		++on_pt;
	}
}
	
double Correlator::getScaleFactor() //returns units of mW
{
	return 2.0*(1000.0*0.5*0.5/(256.0*256.0))*(1.0/(((double)fft_size*fft_size)))/50.0;
}
	
void Waterfall::saveSubWaterfall(int fstart,int fend,int tstart,int tend,double time_offset,string fname)
{
	if(fstart<0) fstart=0;
	if(tstart<0) tstart=0;
	if(fend>=npoints_f) fend=npoints_f-1;
	if(tend>=npoints_t) tend=npoints_t-1;
	cout << fstart << " to " << fend << " and " << tstart << " to " << tend << endl;
	ofstream fout(fname.c_str());
	for(int i=fstart;i<=fend;i++)
	{
	for(int j=tstart;j<=tend;j++)
	{
		fout << i << " " << j << " ";
		fout << ((double)i)*freq_step << " " << ((double)j)*time_step+time_offset << " ";
		double a=data[getIndex(i,j)][0];
		double b=data[getIndex(i,j)][1];
		//fout << data[getIndex(i,j)];
		fout << sqrt(a*a+b*b);
		fout << endl;
	}
	fout << endl;
	}
	fout.close();
}
	
double Waterfall::get_power_squared(int index)
{
	return data[index][0]*data[index][0]+data[index][1]*data[index][1];
}

PowerWaterfall::PowerWaterfall(int npt,int npf,double fstep,double tstep)
{
	data=NULL;
	init(npt,npf,fstep,tstep);
}
	
void PowerWaterfall::init(int npt,int npf,double fstep,double tstep)
{
	npoints_t=npt;
	npoints_f=npf;
	freq_step=fstep;
	time_step=tstep;
	data=new double[npoints_t*npoints_f];
}
	
PowerWaterfall::~PowerWaterfall()
{
	delete data;
}

	
//------------Correlator-----------
	
void Correlator::init(const Monarch *egg,int myfftsize)
{
	fft_size=myfftsize;
    egg->ReadHeader();
    const MonarchHeader *eggheader=egg->GetHeader();
    sampling_rate_mhz=eggheader->GetAcqRate();
    record_size=eggheader->GetRecordSize();
	//decide the optimal size for ffts and allocate memory
    if(record_size<fft_size) {
    	fprintf(stderr,"fft size larger than record size.  make fft size smaller. aborting");
    }
    nffts_per_event=record_size/fft_size;
    fft_output_size=fft_size/2+1;
    fft_input=(float*)fftwf_malloc(sizeof(float)*nffts_per_event*fft_size);
    fft_output_1=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
    fft_output_2=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*fft_output_size*nffts_per_event);
	output_waterfall.init(nffts_per_event,fft_output_size,((double)sampling_rate_mhz*1e6)/((double)fft_size),((double)fft_size)/((double)sampling_rate_mhz*1e6));

    //create the fft plan
    fft_plan_1=fftwf_plan_many_dft_r2c(1,&fft_size,nffts_per_event,fft_input,NULL,1,fft_size,fft_output_1,NULL,1,fft_output_size,FFTW_ESTIMATE);
    fft_plan_2=fftwf_plan_many_dft_r2c(1,&fft_size,nffts_per_event,fft_input,NULL,1,fft_size,fft_output_2,NULL,1,fft_output_size,FFTW_ESTIMATE);


}
	
void Correlator::process_channel1(const Monarch *egg)
{
	const MonarchRecord *event_1=egg->GetRecordOne();
	//convert data to floats
   	for(int i=0;i<record_size;i++)
   		fft_input[i]=(float)(event_1->fDataPtr[i])-128.0;
   	//perform the ffts
   	fftwf_execute(fft_plan_1);
}

void Correlator::process_channel2(const Monarch *egg)
{
	const MonarchRecord *event_2=egg->GetRecordTwo();
	//convert data to floats
	for(int i=0;i<record_size;i++)
   		fft_input[i]=(float)(event_2->fDataPtr[i])-128.0;
   	//perform the ffts
    fftwf_execute(fft_plan_2);
}
	
void Correlator::process_event(const Monarch *egg)
{
	process_channel1(egg);
	process_channel2(egg);
	//make the correlation
	int on_pt=0;
 	for(int i=0;i<nffts_per_event;i++)
   	for(int j=0;j<fft_output_size;j++) {
		output_waterfall.data[on_pt][0]=fft_output_1[on_pt][0]*fft_output_2[on_pt][0]+fft_output_1[on_pt][1]*fft_output_2[on_pt][1];
    	output_waterfall.data[on_pt][1]=fft_output_1[on_pt][0]*fft_output_2[on_pt][1]-fft_output_1[on_pt][1]*fft_output_2[on_pt][0];
		++on_pt;
	}
}

double Correlator::getPowerChannel1(int index)
{
	return fft_output_1[index][0]*fft_output_1[index][0]+
		   fft_output_1[index][1]*fft_output_1[index][1];
}

double Correlator::getPowerChannel2(int index)
{
	return fft_output_2[index][0]*fft_output_2[index][0]+
		   fft_output_2[index][1]*fft_output_2[index][1];
}


#include "Waterfall.hh"
#include "correline_utils.hh"
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

string input_eggname;
string background_filename;
string convolution_filename;
string candidate_list_filename;
int tbins=20;
int fbins=20;

Background background;

void print_usage(); //print out instructions
int handle_options(int argc,char *argv[]); //command line options

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
	if(background_filename=="") {
		cerr << "no background file name given" << endl;
		return -1;
	}
//	if(convolution_filename=="") {
//		cerr << "no convolution file name given" << endl;
//		return -1;
//	}
	if(candidate_list_filename=="") {
		cerr << "no candidate list given " << endl;
		return -1;
	}
	//----load files ----
	if(!background.load(background_filename)) {
		return -1;
	}
	vector<Candidate> candidates=loadCandidateFile(candidate_list_filename);
	if(candidates.size()==0) {
		cerr << "no candidates to view" << endl;
		return -1;
	}
	int fft_size=(background.bg_size-1)*2;
    const Monarch *egg=Monarch::OpenForReading(std::string(input_eggname));
	Correlator correlator;
	correlator.init(egg,fft_size);
	//get candidate events
	for(unsigned int i=0;i<candidates.size();i++) {
		candidates[i].event_no=candidates[i].getEventNo(correlator.sampling_rate_mhz,correlator.record_size);
	}
	
	int on_event=0;
	//Start reading egg
	double total_time=0;
    while(egg->ReadRecord()) {
		//see if there are any candidates in this event
		bool candidate_here=false;
		for(unsigned int i=0;i<candidates.size();i++) {
			if(candidates[i].event_no==on_event)
				candidate_here=true;
		}
		if(candidate_here) {
			correlator.process_event(egg);
			for(unsigned int i=0;i<candidates.size();i++) {
				if(candidates[i].event_no==on_event) {
					cout << "on frequency " << candidates[i].frequency << endl;
					cout << "fstep is " << correlator.output_waterfall.freq_step << endl;

					int fbin=correlator.output_waterfall.getFrequencyBin(candidates[i].frequency*1e6);
					int tbin=correlator.output_waterfall.getTimeBin(candidates[i].time-total_time);
					cout << "fbin " << fbin << " tbin " << tbin << endl;
					stringstream onames;
					onames << "candidate_" << i << ".txt";
					correlator.output_waterfall.saveSubWaterfall(fbin-fbins/2,fbin+fbins/2,tbin-tbins/2,tbin+tbins/2,total_time,onames.str());
				}
			}
		}
		on_event++;
		total_time+=correlator.getEventDuration();
	}
}

void print_usage()
{
	cout << "view_candidate" << endl;
	cout << "prints out a window around a candidate" << endl;
	cout << "options: " << endl;
	cout << "    -i (filename) sets the input egg file MANDATORY" << endl;
	cout << "    -b (filename) sets the ascii background correlation spectrum MANDATORY" << endl;
	cout << "    -l (filename) sets the list of candidates MANDATORY" << endl;
	cout << "    -c (filename) sets the convolution shape" << endl;
	cout << "    -f (integer) sets the number of frequency bins around the candidate to plot" << endl;
	cout << "    -t (integer) sets the number of time bins around the candidate to plot" << endl;
}

int handle_options(int argc,char *argv[])
{
	int c;
    const char *okopts="b:i:l:c:f:t:";
    while((c=getopt(argc,argv,okopts))!=-1)
	switch(c)
	{
		case 'i':
			input_eggname=string(optarg);
			break;
		case 'b':
			background_filename=string(optarg);
			break;
		case 'l':
			candidate_list_filename=string(optarg);
			break;
		case 'c':
			convolution_filename=string(optarg);
			break;
		case 'f':
			fbins=atoi(optarg);
			break;
		case 't':
			tbins=atoi(optarg);
			break;
	}
	return optind;
}

MONARCHLOC := /usr/local
CFLAGS := -Wall -I $(MONARCHLOC)/include -g
LIBS := -L$(MONARCHLOC)/lib -lMonarchCore -lfftw3f -lfftw3f_threads -lpthread -lm -lprotobuf

all: dpph_search powerline correline correline_electron_huntress view_candidate powerline_twochannel powerline_elektronjager

powerline: powerline.c
	g++ $(CFLAGS) -c powerline.c -o powerline.o
	g++ $(CFLAGS) -o powerline $(LIBS) powerline.o

sweepline: sweepline.c
	g++ $(CFLAGS) -c sweepline.c -o sweepline.o
	g++ $(CFLAGS) -o sweepline  $(LIBS) sweepline.o

view_candidate: view_candidate.cc Waterfall.o correline_utils.o
	g++ $(CFLAGS) -c view_candidate.cc -o view_candidate.o
	g++ $(CFLAGS) -o view_candidate  $(LIBS) view_candidate.o Waterfall.o correline_utils.o
	
correline: correline.cc Waterfall.o
	g++ $(CFLAGS) -c correline.cc -o correline.o
	g++ $(CFLAGS) -o correline  $(LIBS) correline.o Waterfall.o

powerline_twochannel: powerline_twochannel.cc Waterfall.o
	g++ $(CFLAGS) -c powerline_twochannel.cc -o powerline_twochannel.o
	g++ $(CFLAGS) -o powerline_twochannel  $(LIBS) powerline_twochannel.o Waterfall.o

dpph_search: dpph_search.cc Waterfall.o
	g++ $(CFLAGS) -c dpph_search.cc -o dpph_search.o
	g++ $(CFLAGS) -o sweepline  $(LIBS) dpph_search.o Waterfall.o

data_simulator: data_simulator.cc correline_utils.o
	g++ $(CFLAGS) -c data_simulator.cc -o data_simulator.o
	g++ $(CFLAGS) -o data_simulator  $(LIBS) data_simulator.o correline_utils.o

powerline_elektronjager: powerline_elektronjager.cc Waterfall.o Histogram.o correline_utils.o
	g++ $(CFLAGS) -c powerline_elektronjager.cc -o powerline_elektronjager.o
	g++ $(CFLAGS) -o powerline_elektronjager  $(LIBS) powerline_elektronjager.o Waterfall.o Histogram.o correline_utils.o


correline_electron_huntress: correline_electron_huntress.cc Waterfall.o Histogram.o correline_utils.o
	g++ $(CFLAGS) -c correline_electron_huntress.cc -o correline_electron_huntress.o
	g++ $(CFLAGS) -o correline_electron_huntress  $(LIBS) correline_electron_huntress.o Waterfall.o Histogram.o correline_utils.o

Waterfall.o: Waterfall.cc Waterfall.hh
	g++ $(CFLAGS) -c Waterfall.cc -o Waterfall.o

Histogram.o: Histogram.cc Histogram.hh
	g++ $(CFLAGS) -c Histogram.cc -o Histogram.o

correline_utils.o: correline_utils.cc correline_utils.hh
	g++ $(CFLAGS) -c correline_utils.cc -o correline_utils.o

clean:
	rm -f *.o
	rm -f powerline
	rm -f sweepline

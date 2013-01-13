MONARCHLOC := /usr/local
CFLAGS := -Wall -I $(MONARCHLOC)/include -g
LIBS := -L$(MONARCHLOC)/lib -lMonarchCore -lfftw3f -lfftw3f_threads -lpthread -lm -lprotobuf

all: dpph_search powerline correline

powerline: powerline.c
	g++ $(CFLAGS) -c powerline.c -o powerline.o
	g++ $(CFLAGS) -o powerline $(LIBS) powerline.o

sweepline: sweepline.c
	g++ $(CFLAGS) -c sweepline.c -o sweepline.o
	g++ $(CFLAGS) -o sweepline  $(LIBS) sweepline.o
	
correline: correline.cc Waterfall.o
	g++ $(CFLAGS) -c correline.cc -o correline.o
	g++ $(CFLAGS) -o correline  $(LIBS) correline.o Waterfall.o

dpph_search: dpph_search.cc Waterfall.o
	g++ $(CFLAGS) -c dpph_search.cc -o dpph_search.o
	g++ $(CFLAGS) -o sweepline  $(LIBS) dpph_search.o Waterfall.o

Waterfall.o: Waterfall.cc Waterfall.hh
	g++ $(CFLAGS) -c Waterfall.cc -o Waterfall.o

#//correline_electron_huntress: correline_electron_huntress.cc Histogram.o
#	g++ $(CFLAGS) -c correline_electron_huntress.cc -o correline_electron_huntress.o
#	g++ $(CFLAGS) -o correline_electron_huntress  $(LIBS) correline_electron_huntress.o Histogram.o

Histogram.o: Histogram.cc Histogram.hh
	g++ $(CFLAGS) -c Histogram.cc -o Histogram.o

clean:
	rm -f powerline
	rm -f sweepline

MONARCHLOC := /usr/local
CFLAGS := -Wall -I $(MONARCHLOC)/include -g
LIBS := -L$(MONARCHLOC)/lib -lmonarch -lfftw3f -lfftw3f_threads -lmxml -lpthread -lm 

all: powerline sweepline

powerline: powerline.c
	g++ $(CFLAGS) -c powerline.c -o powerline.o
	g++ $(CFLAGS) -o powerline $(LIBS) powerline.o

sweepline: sweepline.c
	g++ $(CFLAGS) -c sweepline.c -o sweepline.o
	g++ $(CFLAGS) -o sweepline  $(LIBS) sweepline.o

clean:
	rm -f powerline
	rm -f sweepline

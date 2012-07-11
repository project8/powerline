MONARCHLOC := ../monarch
CFLAGS := -Wall -I $(MONARCHLOC) -g
LIBS := -lfftw3f -lfftw3f_threads -lmxml -lpthread -lm

all: powerline sweepline

powerline: powerline.c
	gcc $(CFLAGS) -o powerline powerline.c $(MONARCHLOC)/monarch.o $(LIBS)

sweepline: sweepline.c
	gcc $(CFLAGS) -o sweepline sweepline.c $(MONARCHLOC)/monarch.o $(LIBS)

clean:
	rm -f powerline
	rm -f sweepline

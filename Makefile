MONARCHLOC := ../monarch
CFLAGS := -Wall -I $(MONARCHLOC) -g
LIBS := -lfftw3f -lfftw3f_threads -lmxml -lpthread

powerline: powerline.c
	gcc $(CFLAGS) -o powerline powerline.c $(MONARCHLOC)/monarch.o $(LIBS)

clean:
	rm -f powerline

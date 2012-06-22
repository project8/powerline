DESCRIPTION:
powerline
Gray Rybka 6/21/2012
Takes a P8 egg file time series and makes a power spectrum out of the first few chunks, printing it to standard out.  Output units are in milliwatts.  Can be in ASCII or binary format, binary is in 64 bit C doubles.

COMPILING:

set location to monarch libraries in Makefile as MONARCHLOC
make

USAGE:

powerline
prints out a power spectrum from an egg fileUsage: powerline [options] (input egg file)
  options:
  -b     sets output to binary (default ASCII)
  -f (integer)  sets the number of points in the fft, default 1024 (about 10 seconds)
  -n (integer)  sets the maximum number of events to scan (default 64)


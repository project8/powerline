#!/usr/bin/env python
import sys
import math
import datetime

def load_histo(filename):
	retval={'x':[],'y':[]}
	thefile=open(filename,"r").readlines()
	for line in thefile:
		elems=line.split(' ')
		retval['x'].append(elems[0])
		retval['y'].append(elems[1])
	return retval

def main():
	sum=None
	for j in range(len(sys.argv)):
		if j==0:
			continue
		arg=sys.argv[j]
		rv=load_histo(arg)
		if sum==None:
			sum=rv
		else:
			for i in range(len(sum['y'])):
				sum['y'][i]=float(sum['y'][i])+float(rv['y'][i])
	for i in range(len(sum['y'])):
		print "%s %g" % (sum['x'][i],float(sum['y'][i]))
#		print sum['x'][i]," ",sum['y'][i]

if __name__=='__main__':
	main()


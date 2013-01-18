#!/usr/bin/env python
import optparse
import math

def load_histo(filename):
	retval={'x':[],'y':[]}
	thefile=open(filename,"r").readlines()
	for line in thefile:
		elems=line.split(' ')
		retval['x'].append(elems[0])
		retval['y'].append(elems[1])
	return retval


def main():
	p=optparse.OptionParser()
	p.add_option('--a','-a',dest="file1",help="first histogram")
	p.add_option('--b','-b',dest="file2",help="second histogram")
	options,arguments=p.parse_args()
	histo1=load_histo(options.file1)
	histo2=load_histo(options.file2)
	logprobsum=0
	for i in range(len(histo1['y'])):
		y1=float(histo1['y'][i])
		y2=float(histo2['y'][i])
		#figure out expected value of lambda of poisson distribution
		m=y1+1 #close enough
		#log likelyhood of y2
		#prob=math.pow(m,y2)*math.exp(-m)/math.factorial(y2)
		logprobsum=logprobsum+y2*math.log(m)-m-math.lgamma(y2+1)
	print -logprobsum/float(len(histo1['y']))


if __name__=='__main__':
	main()



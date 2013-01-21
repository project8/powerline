#!/usr/bin/env python
#compare_histograms
#histogram a is used to determine the probability distribution 
#of the points.  It is assumed histogram a was drawn from poisson statistics
#the likelyhood of obtaining histogram b is then calculated, assuming 
#possion statistics in each bin.  
#this technique is most effective if histogram a is much larger (i.e. has way
#more samples) than histogram b
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

def getsum(histo):
	sum=0.0
	for i in range(len(histo['y'])):
		sum+=float(histo['y'][i])
	return sum

def main():
	p=optparse.OptionParser()
	p.add_option('--a','-a',dest="file1",help="reference histogram (the bigger one)")
	p.add_option('--b','-b',dest="file2",help="histogram to be compared")
	options,arguments=p.parse_args()
	histo1=load_histo(options.file1)
	histo2=load_histo(options.file2)
	#get the normalization of histogram 1 and 2
	histo1sum=getsum(histo1)
	histo2sum=getsum(histo2)
	logprobsum=0
	for i in range(len(histo1['y'])):
		y1=float(histo1['y'][i])
		y2=float(histo2['y'][i])
		#figure out expected value of lambda of poisson distribution
		m=(y1+1.0)*histo2sum/histo1sum #close enough
		#log likelyhood of y2
		#prob=math.pow(m,y2)*math.exp(-m)/math.factorial(y2)
		logprobsum=logprobsum+y2*math.log(m)-m-math.lgamma(y2+1)
	print -logprobsum/float(len(histo1['y']))

if __name__=='__main__':
	main()



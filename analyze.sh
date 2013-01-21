#!/bin/bash
#analyze (file with list of histogram files)

if [ $# -lt 1 ]
then
	echo "usage: analyze (file with list of histogram files)"
	exit 0
fi

filefile=$1
ffbasename=`basename $filefile | sed 's/\(.*\)\..*$/\1/'`
workdir="2013_16_run"
combinedhisto=$workdir"/combinedhisto.txt"
outfile="$workdir""/""$ffbasename""_analysis_result.txt"
#clear out results file
rm -f $outfile

evenfilelist=""
#make a list of only even tens digit offsets
for file in `cat $filefile`
do
	offset=`echo $file | sed 's/.*\(offset[0-9]*\)_.*/\1/'`
	if echo $offset | egrep -q '(2|4|6|8)'
	then
		evenfilelist="$evenfilelist $file"
	fi
done

echo "evenfilelist is $evenfilelist"
python combine_histograms.py $evenfilelist > $combinedhisto
 
for file in `cat $filefile`
do
	#get LO
	lo=`echo $file | sed 's/.*LO\([0-9]*\)_.*/\1/'`
	if [ "$lo" == "$file" ]
	then
		lo=0
	fi
	#get offset
	offset=`echo $file | sed 's/.*offset\([0-9]*\)_.*/\1/'`
	total_offset=$(($lo+$offset))
	echo "lo is $lo"
	num=`python compare_histograms.py -a $combinedhisto -b $file `
	echo "$total_offset $num" >> $outfile
done
sort -n -k 1 $outfile > temp
cp temp $outfile

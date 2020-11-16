#!/bin/sh

mkdir FastScape_$1
cp FastScape FastScape_$1/.
cp *.f90 FastScape_$1/.
rm FastScape_$1/FastScape_NA.f90
rm FastScape_$1/ExtractDistributionFromETOPO1.f90
cp *.f FastScape_$1/.
rm FastScape_$1/na.f
rm FastScape_$1/na_FastScape.f
cp Makefile* FastScape_$1/.
cp FastScape.sh FastScape_$1/.
cp FastScapeView.sh FastScape_$1/.
cp uplift_* FastScape_$1/.
cp initial_topography_* FastScape_$1/.
cp flexure_* FastScape_$1/.
cp precipitation_* FastScape_$1/.
cp tarngo.sh FastScape_$1/.

if [ $# = 2 ]
then
	if [ $2 = GUI ]
	then
	        cp -r FastScape.app FastScape_$1/.
		cp -r pictures FastScape_$1/.
	fi
fi

cp FastScape\ User\ Guide\ Dist.pdf FastScape_$1/.

mkdir FastScape_$1/REF01 FastScape_$1/REF02 FastScape_$1/REF03 FastScape_$1/REF04
mkdir FastScape_$1/REF05 FastScape_$1/REF06 FastScape_$1/REF07
cp REF01/FastScape.in FastScape_$1/REF01/.
cp REF02/FastScape.in FastScape_$1/REF02/.
cp REF03/FastScape.in FastScape_$1/REF03/.
cp REF04/FastScape.in FastScape_$1/REF04/.
cp REF05/FastScape.in FastScape_$1/REF05/.
cp REF06/FastScape.in FastScape_$1/REF06/.
cp REF07/FastScape.in FastScape_$1/REF07/.

tar czvf FastScape_$1.tar.gz FastScape_$1

rm -r FastScape_$1

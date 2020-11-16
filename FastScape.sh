#!/bin/sh -e

EXPERT=0
PREVIEW=0

while getopts ep o
do      case "$o" in
        e)      EXPERT=1;;
        p)      PREVIEW=1;;
        [?])    echo >&2 "Usage: $0 [-ep] DIR_NAME"
                exit 1;;
        esac
done

RUN=$2
if [ $# = 1 ]
then
	RUN=$1
fi

if [ ! -f $RUN/FastScape.in ]
then
	if [ $EXPERT = 1 ]
	then
		cat uplift_start > uplift.f90
		cat uplift_stop >> uplift.f90
		cat initial_topography_start > initial_topography.f90
		cat initial_topography_stop >> initial_topography.f90
		cat flexure_start > flexure.f90
		cat flexure_stop >> flexure.f90
		cat precipitation_start > precipitation.f90
		cat precipitation_stop >> precipitation.f90
		make
	fi
	./FastScape $RUN
fi

if [ $EXPERT = 1 ]
then
  echo "Expert mode"

  mkdir -p $RUN/src
  cp -p *.f $RUN/src
  cp -p *.f90 $RUN/src
  cp Makefile* $RUN/src

  cat uplift_start > $RUN/src/uplift.f90
  awk /uplift_start/,/uplift_stop/ $RUN/FastScape.in >> $RUN/src/uplift.f90
  cat uplift_stop >> $RUN/src/uplift.f90

  cat initial_topography_start > $RUN/src/initial_topography.f90
  awk /initial_topography_start/,/initial_topography_stop/ $RUN/FastScape.in >> $RUN/src/initial_topography.f90
  cat initial_topography_stop >> $RUN/src/initial_topography.f90

  cat flexure_start > $RUN/src/flexure.f90
  awk /dynamic_topography_start/,/dynamic_topography_stop/ $RUN/FastScape.in >> $RUN/src/flexure.f90
  cat flexure_stop >> $RUN/src/flexure.f90

  cat precipitation_start > $RUN/src/precipitation.f90
  awk /precipitation_start/,/precipitation_stop/ $RUN/FastScape.in >> $RUN/src/precipitation.f90
  cat precipitation_stop >> $RUN/src/precipitation.f90

  cd $RUN/src
  make
  cd ../..

fi

RESTART="$(grep restart $RUN/FastScape.in | grep 1 | wc -l)"

if [ $RESTART = 0 ]
then
        rm -f $RUN/*.bmp
        rm -f $RUN/*.pdf
        rm -f $RUN/*.txt
        rm -f $RUN/*.vtk
	rm -f $RUN/*.dem
	rm -f $RUN/*.hdr
        rm -f $RUN/topo*
        rm -f $RUN/temp*
        rm -f $RUN/uplift*
        echo "Previous output erased"
fi

rm -f $RUN/Metric.R

if [ $EXPERT = 0 ]
then
mkdir -p $RUN/src
cp ./FastScape $RUN/src/.
fi

./$RUN/src/FastScape $RUN

if [ -f $RUN/Metric.R ]
then
	R -q -f $RUN/Metric.R > $RUN/R-report.txt
fi

if [ $PREVIEW = 1 ]
then
	echo "Previewing..."
	qlmanage -p $RUN/Topo*.bmp >& /dev/null
	qlmanage -p $RUN/Erosion*.bmp >& /dev/null
	qlmanage -p $RUN/Sediment*.bmp >& /dev/null
	qlmanage -p $RUN/Discharge*.bmp >& /dev/null
	qlmanage -p $RUN/Catchment*.bmp >& /dev/null
	qlmanage -p $RUN/Precip*.bmp >& /dev/null
	qlmanage -p $RUN/Metric*x-av.pdf >& /dev/null
	qlmanage -p $RUN/Metric*y-av.pdf >& /dev/null
	qlmanage -p $RUN/Metric*x-md.pdf >& /dev/null
	qlmanage -p $RUN/Metric*y-md.pdf >& /dev/null
	qlmanage -p $RUN/TopoEvolution.pdf >& /dev/null
	echo "Done"
fi

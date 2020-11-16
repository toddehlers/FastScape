#!/bin/sh -e

qlmanage -p $1/Topo*.bmp >& /dev/null
qlmanage -p $1/Erosion*.bmp >& /dev/null
qlmanage -p $1/Sediment*.bmp >& /dev/null
qlmanage -p $1/Discharge*.bmp >& /dev/null
qlmanage -p $1/Catchment*.bmp >& /dev/null
qlmanage -p $1/Precip*.bmp >& /dev/null
qlmanage -p $1/Metric*x-av.pdf >& /dev/null
qlmanage -p $1/Metric*y-av.pdf >& /dev/null
qlmanage -p $1/TopoEvolution.pdf >& /dev/null

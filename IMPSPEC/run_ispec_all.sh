#!/bin/sh
cd /usr/local/cmod/codes/spectroscopy/vuv/
export shot=$1
./run_ispec_b.sh $1 &
./run_ispec_n.sh $1 &
./run_ispec_o.sh $1 &
./run_ispec_f.sh $1 &
./run_ispec_ne.sh $1 &
./run_ispec_ar.sh $1 &
./run_ispec_mo.sh $1 &

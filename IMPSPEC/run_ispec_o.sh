#!/bin/sh
cd /usr/local/cmod/codes/spectroscopy/vuv/
synchronize_unix loweus_calc_spec $1
synchronize_unix xeus_calc_spec $1
export shot=$1
idl <<EOF
@impspec_ini.bat
shot = long(getenv('shot'))
wait,10
run_impspec,shot,8
exit
EOF

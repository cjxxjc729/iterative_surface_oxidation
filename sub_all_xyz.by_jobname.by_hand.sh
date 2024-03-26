#!/bin/bash
home_dir=$(pwd)
signature=$0
#script_dir=
#tmp_dir=
#mkdir 

#read -p "enter the prefix: " prefix
#read -p "enter the ref pwin file: " f_ref

#--------------------------------------------------------

#n_cif=$(ls -1| grep cif$ | xargs)
job_names=$(ls -1 | grep blank.xyz$ | awk -F '.blank.xyz' '{print $1}' | xargs)


for job_name in $job_names
do
   echo "--- $job_name --------"
   
   ../opt_by_orr_jobname_opt.py $job_name

done

#echo "sleep 2m"
#sleep 2m
#holdon_v1

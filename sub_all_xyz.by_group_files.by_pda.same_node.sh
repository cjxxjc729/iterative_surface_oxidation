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
#job_names=$(ls -1 | grep blank.xyz$ | awk -F '.blank.xyz' '{print $1}' | xargs)

python -c '
import os
def find_and_filter_xyz_files():
    current_dir = "."
    files_and_dirs = os.listdir(current_dir)
    xyz_files = []
    for file in files_and_dirs:
        if file.endswith(".xyz") and not file.endswith("final.xyz"):
            opt_log_file = file[:-4] + "_opt.log"
            if not os.path.exists(opt_log_file):
                xyz_files.append(file)
    with open("fs.t", "w") as f:
        for file in xyz_files:
            f.write(file + "\n")
    print("Found XYZ files:", xyz_files)
if __name__ == "__main__":
    find_and_filter_xyz_files()
'


#ls -1 | grep xyz > fs.t

seperate_list_to_serveral_sublists.sh << EOF
fs.t
20
EOF

group_files=$(ls -1 | grep fs.t.sub_group | xargs)


for group_file in $group_files
do
   echo "--- $group_file --------"
   
   job_num_control 
   hpc_sub_16Gmem "../opt_by_group.py $group_file"

done

echo "sleep 2m"
sleep 2m
holdon_v1

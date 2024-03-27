#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11


#f_dpin='/home/cjx/Desktop/works/my_work/z.Repositories/new_version_test/OC_10M.pb'
#f_dpin='/public1/home/sch0149/script/DPA_suite/OC_300w.pb'
#f_dpin='/public1/home/sch0149/work/subjobs/IMU_leiwang/03.ORR/frozen_model.pb'
#f_ZPE='../ZPE_TS_Eaq/ZPE_TS.json'
#f_Eaq='../ZPE_TS_Eaq/E_aq.json'


import sys
sys.path.append("../")
#sys.path.append("/home/cjx/script/ase_based_constraint_opt_suite/")

#if len(sys.argv) < 2:
#  print("usage: py <prefix>")
#  sys.exit(1)

from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry
import numpy as np
from mimic_functions import *
import os
import time
import re
import glob
import subprocess
import tarfile

from ase.io.lammpsdata import read_lammps_data
from ase.neb import NEB
from ase.optimize import MDMin,BFGS
from deepmd.calculator import DP
from ase.constraints import FixAtoms

from file_operation import *

#input_dic = parse_input_file('input_parameters.json')

#f_ZPE = input_dic['f_ZPE']
#f_Eaq = input_dic['f_Eaq']
#f_reaction_template = input_dic['f_reaction_template']

def collect_traj_4_dpgen(src_directory, dest_directory):
  # 查找以 'traj' 结尾的文件
    traj_files = [os.path.join(root, file)
                  for root, dirs, files in os.walk(src_directory)
                  for file in files if file.endswith("traj")]
    
    print("traj_files=",traj_files)

    if not traj_files:
        return "No 'traj' files found in the source directory."

    # 压缩文件
    max_seq = 0
    existing_files = [file for file in os.listdir(dest_directory) if file.endswith(".tar.gz")]
    for file in existing_files:
        try:
            seq_number = int(file.split('.')[0])
            max_seq = max(max_seq, seq_number)
        except ValueError:
            # 如果文件名不是预期的数字格式，忽略这个文件
            continue
        
    print("max_seq=",max_seq)

    new_file_name = f"{max_seq + 1}.tar.gz"
    new_file_path = os.path.join(dest_directory, new_file_name)

    with tarfile.open(new_file_path, "w:gz") as tar:
        for file in traj_files:
            tar.add(file, arcname=os.path.basename(file))
    
    #return f"Compressed and renamed to {new_file_path}"



def sub_opt_by_jobs():
  
  start_time = time.time()
  input_dic = parse_input_file('../input_parameters.json')
  sub_script = input_dic['opt_script'] 

  # 命令和需要传递给命令的输入
  #command = ['../sub_all_cif.by_jobname.by_pda.same_node.sh']
  #command = ['../sub_all_xyz.by_jobname.by_pda.same_node.sh']
  #command = ['../sub_all_xyz.by_jobname.by_hand.sh']
  command = ['../'+sub_script]
  print("shell command =", command)

  # 执行命令，并通过input参数传递输入
  with open('sub_all_xyz.by_jobname.by_pda.same_node.out','w') as output_file:
    subprocess.run(command, text=True, stdout=output_file, stderr=subprocess.PIPE)

  #print("sleep 120s")
  #time.sleep(120)

  #command = ['holdon_v1']
  #print("shell command =", command)

  # 执行命令，并通过input参数传递输入
  #with open('sub_all_cif.by_jobname.by_pda.same_node.out','a') as output_file:
  #  subprocess.run(command, text=True, stdout=output_file, stderr=subprocess.PIPE)

  end_time = time.time()

  
  save_str_4_model_dev = input_dic['save_str_4_model_dev'] 
  
  if save_str_4_model_dev == 'True':  #如果save_str_4_model_dev == 'True'， 则收集traj文件4后续的DPGEN任

    print("compress traj")
    mkdir_if_not_exists('../traj_coll')
    collect_traj_4_dpgen('./','../traj_coll')
    


  elapsed_time = end_time - start_time  # Calculate the elapsed time

  print("end sub_opt_by_jobs()")
  print("Elapsed time: {:.2f} seconds".format(elapsed_time))

  


if __name__ == "__main__":

  input_dic = parse_input_file('../input_parameters.json')

  f_dpin = input_dic['f_dpin']
  f_ZPE = input_dic['f_ZPE']
  f_Eaq = input_dic['f_Eaq']
  f_reaction_template = input_dic['f_reaction_template']
  max_opt_step = float(input_dic['max_opt_step'])

  jobname=sys.argv[1]
  # Pattern to match files starting with 'C1027.Rblank' and ending with '.cif'
  pattern = './' + jobname +'*.xyz'

  # Using glob.glob() to find all files matching the pattern
  fs_str = glob.glob(pattern)
  print("fs_str=",fs_str)


  for f_str in fs_str:
    print("------------------",f_str,"--------------------")

    atoms = read(f_str)
    prefix = f_str.split('.xyz')[0]

    if os.path.exists(prefix+'.fix_list'):
      print("find fix list ")
      fix_list=np.loadtxt(prefix+'.fix_list',ndmin=1)-1
      c = FixAtoms(indices=fix_list)
      atoms.set_constraint(c)

    atoms.calc = DP(model=f_dpin)
    opt_ini = BFGS(atoms, trajectory=prefix+'_opt.traj',  logfile=prefix+'_opt.log')
    opt_ini.run(fmax=0.05,steps=max_opt_step)

    atoms_final=read(prefix+'_opt.traj')
    write(prefix+'_final.xyz',atoms_final)


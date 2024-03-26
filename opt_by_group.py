#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11

#f_dpin='/public1/home/sch0149/script/DPA_suite/OC_300w.pb'
f_dpin='/public1/home/sch0149/work/subjobs/IMU_leiwang/03.ORR/frozen_model.pb'
f_ZPE='../ZPE_TS_Eaq/ZPE_TS.json'
f_Eaq='../ZPE_TS_Eaq/E_aq.json'


import sys
sys.path.append("/public1/home/sch0149/script/mimic_suite/")
sys.path.append("/public1/home/sch0149/script/ase_based_constraint_opt_suite/")

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

from ase.io.lammpsdata import read_lammps_data
from ase.neb import NEB
from ase.optimize import MDMin,BFGS
from deepmd.calculator import DP
from ase.constraints import FixAtoms


def sub_opt_by_jobs():
  
  start_time = time.time()
 
  # 命令和需要传递给命令的输入
  #command = ['../sub_all_cif.by_jobname.by_pda.same_node.sh']
  command = ['../sub_all_xyz.by_group_files.by_pda.same_node.sh']
  print("shell command =", command)

  # 执行命令，并通过input参数传递输入
  with open('sub_all_xyz.by_group_files.by_pda.same_node.out','w') as output_file:
    subprocess.run(command, text=True, stdout=output_file, stderr=subprocess.PIPE)

  #print("sleep 120s")
  #time.sleep(120)

  #command = ['holdon_v1']
  #print("shell command =", command)

  # 执行命令，并通过input参数传递输入
  #with open('sub_all_cif.by_jobname.by_pda.same_node.out','a') as output_file:
  #  subprocess.run(command, text=True, stdout=output_file, stderr=subprocess.PIPE)

  end_time = time.time()

  elapsed_time = end_time - start_time  # Calculate the elapsed time

  print("end sub_opt_by_jobs()")
  print("Elapsed time: {:.2f} seconds".format(elapsed_time))

  


if __name__ == "__main__":

  #jobname=sys.argv[1]
  # Pattern to match files starting with 'C1027.Rblank' and ending with '.cif'
  #pattern = './' + jobname +'*.xyz'

  group_file = sys.argv[1]

  # Using glob.glob() to find all files matching the pattern
  #fs_str = glob.glob(pattern)
  #print("fs_str=",fs_str)

  fs_str = np.loadtxt(group_file,dtype=str)

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
    opt_ini.run(fmax=0.05,steps=200)

    atoms_final=read(prefix+'_opt.traj')
    write(prefix+'_final.xyz',atoms_final)


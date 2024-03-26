#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11


#f_ZPE='../ZPE_TS_Eaq/ZPE_TS.json'
#f_Eaq='../ZPE_TS_Eaq/E_aq.json'

import json

# Read the JSON data from the file
#with open(f_ZPE, 'r') as file:
#    data_ZPE = json.load(file)

#with open(f_Eaq, 'r') as file:
#    data_Eaq = json.load(file)

import sys
sys.path.append("../")
#sys.path.append("/home/cjx/script/ase_based_constraint_opt_suite/")

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

from gen_FED import *



def cal_dE(FED):

  if len(FED) == 2:
    print("only one step")
    dE=FED[1]-FED[0]
    dE=dE/2
    TDI = 0
    TDTSv = 0
    return dE, TDI, TDTSv    
 
  else:

    #1 make G05
  
    dG_05=[]
    #print(len(FED))
    for k in range(len(FED)-1):
      dG_05.append((FED[k]+FED[k+1])/2)

    #-------- main  pro -------------------------------!
    M_dE_minus = []
    for i in range(len(FED)-1):
      for j in range(i+1,len(FED)-1):
        M_dE_minus = np.append(M_dE_minus,[dG_05[i]-FED[j],i,j])

    M_dE_minus = M_dE_minus.reshape(-1,3)
    max_idx_minus = np.argmax(M_dE_minus[:,0])
    dE_minus = M_dE_minus[max_idx_minus,0]

    M_dE_plus = []
    for j in range(len(FED)-1):
      for i in range(j,len(FED)-1):
        M_dE_plus = np.append(M_dE_plus,[dG_05[i]-FED[j],i,j])

    M_dE_plus = M_dE_plus.reshape(-1,3)
    max_idx_plus = np.argmax(M_dE_plus[:,0])
    dE_plus = M_dE_plus[max_idx_plus,0]

    Gtot=FED[-1]-FED[0]
    dE=max(dE_minus+Gtot,dE_plus)

    if dE_minus+Gtot > dE_plus:
      TDTSv = M_dE_minus[max_idx_minus,1]
      TDI   = M_dE_minus[max_idx_minus,2]
    else:
      TDTSv = M_dE_plus[max_idx_plus,1]
      TDI   = M_dE_plus[max_idx_plus,2]

    TDI   = int(TDI)
    TDTSv = int(TDTSv)

    return dE, TDI, TDTSv


def get_TDI(jobname ,U_vs_RHE, data_ZPE, data_Eaq):

  # Pattern to match files starting with 'C1027.Rblank' and ending with '.cif'
  pattern = './' + jobname +'*.log'

  # Using glob.glob() to find all files matching the pattern
  fs_log = glob.glob(pattern)
  #print("fs_str=",fs_log)

  Gibbs_OER=[]
  GH2O = data_Eaq['H2O'][0]
  GH   = data_Eaq['H+'][0]
  GO2  = data_Eaq['O2'][0]

  dic_Gadb={}

  for f_log in fs_log:

    adb_name = f_log.split('.')[-2].split('_')[0]
    print("f_log = ", f_log)
    with open(f_log,'r') as f:
      E = f.readlines()[-1].split()[-2]
    E = float(E)

    G = E + data_ZPE[adb_name][0] - data_ZPE[adb_name][1] + data_ZPE[adb_name][2]

    dic_Gadb[adb_name] = G

  #print("dic_Gadb=",dic_Gadb)
  f_cat_rea_adbs='../ref_str/cat_rea_adbs.txt'

  FED, adb_list = make_FED_test(dic_Gadb, data_Eaq, f_cat_rea_adbs, U_vs_RHE)


  dEv, TDI_idx, TDTSv_idx = cal_dE(FED)
  TDI = adb_list[TDI_idx]

  return TDI, dEv, FED


def make_FED_test(dic_Gadb, data_Eaq, f_cat_rea_adbs, U):

  #1
  adb_list = np.loadtxt(f_cat_rea_adbs,dtype=str)
  
  Ge= -U

  GH2O = data_Eaq['H2O'][0]
  GH   = data_Eaq['H+'][0]
  GO2  = data_Eaq['O2'][0]
  
  print("dic_Gadb",dic_Gadb)

  Gibbs_OER = []
  Gibbs_OER.append(2*GH2O+dic_Gadb['blank'])
  Gibbs_OER.append(GH2O+dic_Gadb['OH']+GH+Ge)
  Gibbs_OER.append(GH2O+dic_Gadb['O']+2*GH+2*Ge)
  Gibbs_OER.append(dic_Gadb['OOH']+3*GH+3*Ge)
  Gibbs_OER.append(dic_Gadb['blank']+4*GH+4*Ge+GO2)

  Gibbs_OER = np.array(Gibbs_OER)
  FED = Gibbs_OER - Gibbs_OER[0]

  return FED, adb_list



def select_opt_str():

  #------------------
  input_dic = parse_input_file('../input.parameters')

  f_ZPE = input_dic['f_ZPE']
  f_Eaq = input_dic['f_Eaq']
  URHE = input_dic['URHE']
  catalytic_reaction_name=input_dic['catalytic_reaction_name']
  cal_tot_signal = input_dic['cal_tot'] 

  with open(f_ZPE, 'r') as file:
    dic_ZPE = json.load(file)
  with open(f_Eaq, 'r') as file:
    dic_Eaq = json.load(file)

  #--------------

  files = os.listdir('.')

  # Filter out files where z='blank' and s='cif'
  #filtered_files = [file for file in files if file.endswith("blank.cif")]
  filtered_files = [file for file in files if file.endswith("blank.xyz")]

  # Extract the x and y values from the filtered files
  xy_values = [file.split('.')[:2] for file in filtered_files]

  # Convert list of lists into a list of tuples to make them hashable for set operations
  xy_tuples = [tuple(xy) for xy in xy_values]

  # Find unique x and unique y values
  unique_x = set(x for x, _ in xy_tuples)
  unique_y = set(y for _, y in xy_tuples)

  if cal_tot_signal == 'False':
    atoms_new_TDI = read('shell.xyz')
  else:
    atoms_new_TDI = Atoms()

  for cluster in unique_x:
    print("----cluster", cluster,"---------------")
    results_role_TDI_dEv = []
    for role in unique_y:
      #if os.path.exists(cluster+'.'+role+'.blank.cif'):
      if os.path.exists(cluster+'.'+role+'.blank.xyz'):
        
        jobname = cluster+'.'+role
        dic_Gadb = get_dic_Gadb(jobname, dic_ZPE)
        FED, adb_list = gen_FED(jobname, catalytic_reaction_name, dic_Gadb, dic_Eaq, URHE) 
        dEv, TDI_idx, TDTSv_idx = cal_dE(FED)
        TDI = adb_list[TDI_idx]    

        
        #TDI, dEv, FED  = get_TDI(f_jobname ,U_vs_RHE, data_ZPE, data_Eaq)
 
        results_role_TDI_dEv = np.append(results_role_TDI_dEv, [role, TDI, dEv])
  
    results_role_TDI_dEv = results_role_TDI_dEv.reshape(-1,3)

    idx_of_reaction_that_actually_happen = np.argmin(results_role_TDI_dEv[:,0])  #dEv更小的那个会发生
    
    i = idx_of_reaction_that_actually_happen
    role = results_role_TDI_dEv[i, 0]
    TDI  = results_role_TDI_dEv[i, 1]
    dEv  = results_role_TDI_dEv[i, 2]

    f_final_xyz = cluster+'.'+role+'.'+TDI+'_final.xyz'
    print("dEv = ", dEv)
    print(f_final_xyz)
    atom_append = read(f_final_xyz)   
    #geometry.get_duplicate_atoms(atom_append, 0.1, delete=True) 

    np.savetxt(cluster+'.'+role+'.FED',FED  ,fmt='%.4f')
    np.savetxt(cluster+'.'+role+'.dEv',[dEv],fmt='%s') 
    np.savetxt(cluster+'.'+role+'.TDI',[TDI],fmt='%s')

    atoms_new_TDI = atoms_new_TDI + atom_append
  
  geometry.get_duplicate_atoms(atoms_new_TDI, 0.4, delete=True)
  atoms_new_TDI_dot = Atoms(atoms_new_TDI.get_chemical_symbols(),positions=atoms_new_TDI.get_positions(),cell=atoms_new_TDI.get_cell())

  return atoms_new_TDI_dot

if __name__ == "__main__":

  atoms_new_TDI = select_opt_str() 

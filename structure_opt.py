#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11

import sys
#sys.path.append("/public1/home/sch0149/script/mimic_suite/")
#sys.path.append("/public1/home/sch0149/ase_based_constraint_opt_suite/")
#sys.path.append("/public1/home/sch0149/script/mmpp/")

from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry
import numpy as np
from mimic_functions import *
from file_operation import *

import os
import time
import re
import shutil

def make_cfid(center_atom_id):
  a=os.path.exists("c"+str(center_atom_id))
  if a==False:
    b=os.getcwd()
    os.mkdir(b+"/c"+str(center_atom_id))


def put_atom_to_cell_center_with_large_vac(atoms):

  maxx=max(atoms.get_positions()[:,0])
  minx=min(atoms.get_positions()[:,0])
  maxy=max(atoms.get_positions()[:,1])
  miny=min(atoms.get_positions()[:,1])
  maxz=max(atoms.get_positions()[:,2])
  minz=min(atoms.get_positions()[:,2])
  lenx=maxx-minx
  leny=maxy-miny
  lenz=maxz-minz
  length=max(3*lenx,3*leny,3*lenz,50)

  #pos=atoms.get_positions()+[length/2,length/2,length/2]
  atoms.set_cell([(length,0,0),(0,length,0),(0,0,length)])
  #atoms.center()

def up_to_z(atoms_nl):

  print("calculae the max dis in ref atoms\n")

  D,Dlen=geometry.get_distances(atoms_nl.get_positions())
  max_dis_in_ref=np.max(Dlen)

  print("set cutoff. cutoff is set based on the size of the ref fragment. it is the largest distance in ref atoms plus 0.1")
  print("the largest distance is "+str(max_dis_in_ref))

  cutoff_for_box=max_dis_in_ref+0.1
  side_length_of_box=cutoff_for_box+2


  c_atm_id=0
  print("\n-----------------------shift the cor stick up to the z axis ----------------------------\n")
  atoms_nl_rotted,M_rot=rot_to_z_axis_surface_version(atoms_nl, c_atm_id, side_length_of_box)
  atoms_nl_rotted.center()
  print("M_rot")
  print(M_rot)

  print("\n---------------end of shift the cor stick up to the z axis---------------------\n")


  return atoms_nl_rotted, M_rot


def gen_shell_reover_atoms(atoms,c_atm_id, nl,atoms_nl_rotted, M_rot):


  print("\n-------------------------- we here generate a shell atom----------------------------\n")

  v1=-atoms.get_positions()[c_atm_id]
  v2=atoms_nl_rotted.get_positions()[0]
  cors=atoms.get_positions()
  rotted_cors=(cors+v1)*M_rot+v2
  atoms.set_positions(rotted_cors)


  indices, offsets = nl.get_neighbors(c_atm_id)
  indices=np.array(indices)

  full_indices=list(range(len(atoms)))
  full_indices=np.array(full_indices)

  mask=np.isin(full_indices,indices)

  shell_index = full_indices[~mask]
  shell_index = list(shell_index)


  if len(shell_index) > 0:
    atoms_shell=atoms[shell_index]
  else:
    print("no shell valuem, no recover.xyz")


def gen_fix_list(atoms_nl, cutoff_i):

  for atm_id in range(len(atoms_nl)):
    if atoms_nl.get_distance(0,atm_id) >= cutoff_i:
      break
 
  fix_list=[x for x in range(len(atoms_nl)) if x >= atm_id]
  fix_list=np.array(fix_list)+1

  return fix_list



def opt_str_ini(f_str, c_atm_ids):

  cutoff_o=7
  cutoff_i=4
  
  atoms=read(f_str)    
  
  atoms_nls,nl=gen_atoms_nl_4_atm_ids(atoms, c_atm_ids, set_cutoff=cutoff_o)
  
  dir_name='str_opt_center'
  mkdir(dir_name)
  
  atm_ids_sel_tot=[]
  for c_atm_id in c_atm_ids:
  
    #1
    atoms_core = atoms_nls[c_atm_id][1:]
     
    #2
    fix_list = gen_fix_list(atoms_core, cutoff_i)
    np.savetxt(dir_name+'/'+str(c_atm_id+1)+'_core.fix_list',fix_list,fmt='%d')
  
    #3
    atm_ids_core, offset = nl.get_neighbors(c_atm_id)
    write(dir_name+'/'+str(c_atm_id+1)+'_core.cif',atoms_core)
  
    atm_ids_sel_tot=np.append(atm_ids_sel_tot,atm_ids_core)
  
  
  atm_ids_shell = [x for x in range(len(atoms)) if x not in atm_ids_sel_tot]
  atoms_shell = atoms[atm_ids_shell]
  write(dir_name+'/'+'atoms_shell.cif',atoms_shell)



if __name__ == "__main__":
   
  cutoff_o=7
  cutoff_i=4

  atoms=read('1.xyz')    
  c_atm_ids=[990,3961]

  atoms_nls,nl=gen_atoms_nl_4_atm_ids(atoms, c_atm_ids, set_cutoff=cutoff_o)

  dir_name='str_opt_center'
  mkdir(dir_name)

  atm_ids_sel_tot=[]
  for c_atm_id in c_atm_ids:
  
    #1
    atoms_core = atoms_nls[c_atm_id][1:]
     
    #2
    fix_list = gen_fix_list(atoms_core, cutoff_i)
    np.savetxt(dir_name+'/'+str(c_atm_id+1)+'_core.fix_list',fix_list,fmt='%d')

    #3
    atm_ids_core, offset = nl.get_neighbors(c_atm_id)
    write(dir_name+'/'+str(c_atm_id+1)+'_core.cif',atoms_core)

    atm_ids_sel_tot=np.append(atm_ids_sel_tot,atm_ids_core)


  atm_ids_shell = [x for x in range(len(atoms)) if x not in atm_ids_sel_tot]
  atoms_shell = atoms[atm_ids_shell]
  write(dir_name+'/'+'atoms_shell.cif',atoms_shell)
  
    


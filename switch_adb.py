#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11


import sys
#s#ys.path.append("/public1/home/sch0149/script/mimic_suite/")
#sys.path.append("/public1/home/sch0149/ase_based_constraint_opt_suite/")
#sys.path.append("/public1/home/sch0149/script/mmpp/")


from scipy.spatial.transform import Rotation as R
from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry
import numpy as np
from mimic_functions import *
#from mmpp_functions import *
#from project_to_grid_functions import *
import os
import time
import re



#def make_adb_by_switch(atoms, adb_name_from, adb_name_to, dir_4_ref_str):
  
#  '''
#  '''
#  #
#  natm     = len(atoms)
#  natm_adb = len(adb_name_from)
 
#  atoms_adb  = atoms[-natm_adb:]  #atoms of adb
#  atoms_base = atoms[:-natm_adb]  #atoms of base
  
#  f_refstr_adb_from = dir_4_ref_str+'/'+adb_name_from +'.cif'
#  f_refstr_adb_to   = dir_4_ref_str+'/'+adb_name_to   +'.cif'

#  atoms_ref_str_of_adb_from = read(f_refstr_adb_from)
#  atoms_ref_str_of_adb_to   = read(f_refstr_adb_to)

#  #2
#  ref_ad_site_from, ref_ad_atm_from, ref_atoms_adb_from  = read_site_and_keyadb(atoms_ref_str_of_adb_from)
#  ref_ad_site_to,   ref_ad_atm_to  , ref_atoms_adb_to    = read_site_and_keyadb(atoms_ref_str_of_adb_to)  #
  
#  norm_vec = vect_for_most_sparse_locally_v2(atoms, natm-natm_adb)  #通过ad_atm_from元素坐标计算norm_vec, its atm_id is natm-natm_adb
 
#  #通过norm_vec, Ar元素位置（即吸附中心位置）， 将新的吸附质的吸附中心和Ar元素位置重合

   
 

def switch_adb(adb0_name,adb_name,atoms_adb0,atm_ids_adb_of_atoms_adb0,dir_4_ref_str):

  '''
  switch_adb主程序, adb0
  adb0表示原有adb，adb表示要转换的adb。ref_str_of_adb是structur file of吸附质，通常是单原子带一个吸附质的结构。
  如OH.cif这类。这类cif的特点是最后几个原子序号为吸附质，吸附质之上的那个是吸附位，如：
  Fe 0.5 0.5 0.5  #ad_site
  O  0.5 0.5 0.6  #ad_atm
  H  0.5 0.5 0.65
  '''

  #1 ini  
  f_refstr_adb0 = dir_4_ref_str+'/'+adb0_name+'.cif'
  f_refstr_adb  = dir_4_ref_str+'/'+adb_name +'.cif'

  atoms_ref_str_of_adb0 = read(f_refstr_adb0) 
  atoms_ref_str_of_adb  = read(f_refstr_adb) 

  #2
  ref_ad_site0, ref_ad_atm0, ref_atoms_adb0  = read_site_and_keyadb(atoms_ref_str_of_adb0)
  ref_ad_site,  ref_ad_atm,  ref_atoms_adb   = read_site_and_keyadb(atoms_ref_str_of_adb)
  #print("ref_ad_site,  ref_ad_atm,  ref_atoms_adb", ref_ad_site,  ref_ad_atm,  ref_atoms_adb) 
  
  #3
  atoms=atoms_adb0.copy()  #创造一个副本
  cor_ad_atm0 = atoms[atm_ids_adb_of_atoms_adb0[0]].position  #找到吸附位中心ad_atm的坐标
  #print("atm_ids_adb_of_atoms_adb0[0]",atm_ids_adb_of_atoms_adb0[0])
  #print("cor_ad_atm0=",cor_ad_atm0)

  del atoms[atm_ids_adb_of_atoms_adb0]
  print("ref_ad_atm=",ref_ad_atm)

  #if len(ref_ad_atm) != 0:  #case of not blank
     
  atoms = atoms+Atoms('Ar',positions=[cor_ad_atm0])  #删去吸附质并在吸附中心加上Ar元素
  
  norm_vec = vect_for_most_sparse_locally_v2(atoms,len(atoms)-1)  #通过Ar元素坐标计算norm_vec, 固而 atm_id = len(atoms)-1
  
  #通过norm_vec, Ar元素位置（即吸附中心位置）， 将新的吸附质的吸附中心和Ar元素位置重合，将其ref_vec （即[0,0,1]）平移到当前的norm_vec上。
  
  ref_atoms_adb.translate(-ref_atoms_adb[0].position)  #step1:平移至原点
  a,v = calculate_rot_av_from_refvec1_to_refvec2([0,0,1],norm_vec)
  ref_atoms_adb.rotate(a,v)
  ref_atoms_adb.translate(cor_ad_atm0)  #step3:平移至Ar cor
   
  #del atoms_adb0[atm_ids_adb_of_atoms_adb0]
  base_idx = list(range(len(atoms_adb0)))
  base_idx = [x for x in base_idx if x not in atm_ids_adb_of_atoms_adb0]
  atoms_adb = atoms_adb0[base_idx]+ref_atoms_adb

  #else: #blank

  #atoms_adb = atoms

  return atoms_adb



def read_site_and_keyadb(atoms_ref_str_of_adb):
   
  atoms = atoms_ref_str_of_adb  #必须是ref_str文件夹里面的H.cif OH.cif之类的文件

  #1
  eles = atoms.get_chemical_symbols()
  eles = np.array(eles)
  Fe_index = np.where(eles == 'Fe')[0][0]
  adb_index = list(range(int(Fe_index)+1,len(atoms)))
#np.atleast_2d(list(range(int(Fe_index)+1,len(atoms)))) # 使用np.atleast_2d来确保数组至少有两个维度,使得下面可以用positions

  #2
  cor_Fe   = atoms[Fe_index].position
  if len(adb_index) != 0:
    cors_adb = atoms[adb_index].positions
    D,Dlen=geometry.get_distances(cors_adb,cor_Fe)
    atm_id_keyadb = adb_index[np.argmin(Dlen)]

    #print("atm_id_keyadb = ",atm_id_keyadb)
    ad_site = Fe_index
    ad_atm = atm_id_keyadb

    #3
    atoms_adb = atoms[adb_index]  #由于上文已经使用了np.atleast_2d这边出来的一定是Atoms
    #atoms_adb.set_pbc((0, 0, 0))  #cancel pbc
    atoms_adb.cell = np.zeros((3, 3))
    #print("atoms_adb = ", atoms_adb)


  else:  # for adb = blank
    ad_site = Fe_index
    ad_atm = [] 
    atoms_adb = []

  return ad_site, ad_atm, atoms_adb


def calculate_rot_av_from_refvec1_to_refvec2(refvec1,refvec2):
  
  '''
  a=angel
  v=axis
  '''

  A1=refvec1
  A2=refvec2

  # 计算旋转轴（A1和A2的叉积）
  rotation_axis = np.cross(A1, A2) 
  
  # 计算旋转角度
  angle_cos = np.dot(A1, A2) / (np.linalg.norm(A1) * np.linalg.norm(A2))
  angle_rad = np.arccos(angle_cos)
  # 将弧度转换为角度
  angle_deg = np.degrees(angle_rad)

  # 构建四元数表示的旋转
  #rotation = R.from_rotvec(rotation_axis * angle_rad / np.linalg.norm(rotation_axis))

  # 应用旋转
  #rotated_point = rotation.apply(point)

  #cor2 = rotated_point
 
  a = angle_deg
  v = rotation_axis

  return a, v
 
def add_adb(adb1_name,atoms,charact_atm_ids_adb,dir_4_ref_str):
 
  '''
  add_adb主程序, 其为atoms0为blank.cif时的特例
  adb0表示原有adb，adb表示要转换的adb。ref_str_of_adb是structur file of吸附质，通常是单原子带一个吸附质的结构。
  如OH.cif这类。这类cif的特点是最后几个原子序号为吸附质，吸附质之上的那个是吸附位，如：
  Fe 0.5 0.5 0.5  #ad_site
  O  0.5 0.5 0.6  #ad_atm
  H  0.5 0.5 0.65
  '''

  #1 ini  
  f_refstr_adb  = dir_4_ref_str+'/'+adb1_name +'.cif'

  atoms_ref_str_of_adb  = read(f_refstr_adb)
  
  #2
  ref_ad_site,  ref_ad_atm,  ref_atoms_adb   = read_site_and_keyadb(atoms_ref_str_of_adb) 
 
  #3

  #3.1 swift c_atm_id to the end
  atoms_adb = atoms.copy()
  #del atoms_adb[charact_atm_ids_adb]
  #atoms_site = atoms[charact_atm_ids_adb]
  #atoms_adb = atoms_adb + atoms_site

  #3.2 add adb, only work if not blank
  if len(ref_atoms_adb) != 0:  #case of not blank

    c_atm_id = charact_atm_ids_adb[0]
    #c_atm_id = len(atoms_adb)-1 # since c_atm has been swifted to the end, c_atm_id should equal len(atoms)-1

    norm_vec = vect_for_most_sparse_locally_v2(atoms_adb,c_atm_id)  #通过c_atm_id元素坐标计算norm_vec
    
    ref_atoms_adb = atoms_ref_str_of_adb[[ref_ad_site]]+ ref_atoms_adb
    ref_atoms_adb.cell = np.zeros((3,3))    #更改一下ref_atom_adb 使之囊括吸附位

    #通过norm_vec, 即吸附中心位置， 将新的吸附质的吸附中心和Ar元素位置重合，将其ref_vec （即[0,0,1]）平移到当前的norm_vec上。  
    ref_atoms_adb.translate(-ref_atoms_adb[0].position)  #step1:平移至原点
    a,v = calculate_rot_av_from_refvec1_to_refvec2([0,0,1],norm_vec)
    ref_atoms_adb.rotate(a,v)
    ref_atoms_adb.translate(atoms_adb[c_atm_id].position)  #step3:平移至c_atm_id
   
    del ref_atoms_adb[0]
    atoms_adb = atoms_adb+ref_atoms_adb

  return atoms_adb


def get_c_atm_id_from_adb(atoms, charact_atm_ids_adb, nl):

  indices, offset = nl.get_neighbors(charact_atm_ids_adb[0]) 
  atm_ids = sort_the_indices(atoms, charact_atm_ids_adb[0], indices)[:,0]
  atm_ids = [int(x) for x in atm_ids]
  for atm_id in atm_ids:
    if atm_id not in charact_atm_ids_adb:
      c_atm_id = atm_id
      break
  return c_atm_id

if __name__ == "__main__":
  
  #input
  adb0_name='OH'
  adb1_name='blank'

  atoms_adb0=read('1.xyz')    
  charact_atm_ids_adb=[125,114]
  nl = ini_nl(atoms_adb0,set_cutoff=4.5)
  
  dir_4_ref_str='/home/cjx/Desktop/test/surf_mesh_test/ORR_test/ref_str'
 
  #output
  if adb0_name != 'blank':    
    c_atm_id = get_c_atm_id_from_adb(atoms_adb0, charact_atm_ids_adb, nl)
    #atoms_adb1 = switch_adb(adb0_name,adb1_name,atoms_adb0,charact_atm_ids_adb,dir_4_ref_str)
  else:
    c_atm_id = charact_atm_ids_adb
  
  atoms_adb1 = add_adb(adb1_name,atoms_adb0,c_atm_id,dir_4_ref_str)

  write('atoms_adb1.xyz',atoms_adb1)
  


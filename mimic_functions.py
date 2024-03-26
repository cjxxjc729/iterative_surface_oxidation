#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11

import sys
import numpy as np
from ase.io import write
from ase.io import read
#sys.path.append("/home/cjx/script/ase_based_constraint_opt_suite/")
#from ase.geometry.analysis import Analysis
from ase.geometry import get_duplicate_atoms,get_distances,get_angles
from ase import neighborlist
#from project_to_grid_functions import *
from ase import Atoms
import os
import random
import time
import shutil

def ini_nl(atoms,set_cutoff):

  print("entering ini_nl")  

  start_time = time.time()

  cutOff=[]
  for i in range(len(atoms)):
    cutOff.append(set_cutoff/2)
  cutOff=np.array(cutOff)
  nl = neighborlist.NeighborList(cutOff, skin=0, bothways=True)
  nl.update(atoms)

  end_time = time.time()

  elapsed_time = end_time - start_time  # Calculate the elapsed time

  print("end ini_nl")
  print("Elapsed time: {:.2f} seconds".format(elapsed_time))

  return nl


def get_nei_indices_v2(atoms,atm_id,set_cutoff,nl):
  
  indices, offsets = nl.get_neighbors(atm_id)
  atm_nl=sort_the_indices(atoms,atm_id,indices)
  indices=atm_nl[:,0]
  indices=indices.astype(int)  #this will bring in a sort to the neibor list
  indices=list(indices)

  return indices


def get_nei_indices(atoms,atm_id,set_cutoff):
  cutOff = neighborlist.natural_cutoffs(atoms)
  for i in range(len(cutOff)):
    cutOff[i]=set_cutoff/2
  nl = neighborlist.NeighborList(cutOff, skin=0, bothways=True)
  nl.update(atoms)
  indices, offsets = nl.get_neighbors(atm_id)
  atm_nl=sort_the_indices(atoms,atm_id,indices)
  indices=atm_nl[:,0]
  indices=indices.astype(int)  #this will bring in a sort to the neibor list
  indices=list(indices)
  return indices

def get_nei_indices_from_nature_bond(atoms,atm_id):
  cutOff = neighborlist.natural_cutoffs(atoms)
  nl = neighborlist.NeighborList(cutOff, skin=0.3, bothways=True)
  nl.update(atoms)
  indices, offsets = nl.get_neighbors(atm_id)
  atm_nl=sort_the_indices(atoms,atm_id,indices)
  indices=atm_nl[:,0]
  indices=indices.astype(int)  #this will bring in a sort to the neibor list
  indices=list(indices)
  return indices

def get_neishell_indices(atoms,atm_id,dis_range):
  cuto=dis_range[1]
  cuti=dis_range[0]
  nei_ids_o=get_nei_indices(atoms,atm_id,cuto)
  nei_ids_i=get_nei_indices(atoms,atm_id,cuti)
  for atm_id in nei_ids_i:
    nei_ids_o.remove(atm_id)
  nei_ids_io=nei_ids_o  
  return nei_ids_io

def count_rdf_for_givenid(atoms,atm_id,dis_cutoff=10,ds=0.5):
  n_ds=int(dis_cutoff/ds)
  rdf=[]
  for i in range(n_ds):
    dis_range=[0.1+i*ds,0.1+i*ds+ds]
    x=np.mean(dis_range)
    print("-------------range of "+str(dis_range)+" ---------------------")
    compose_array=count_neishell_compose_for_givenid(atoms,atm_id,dis_range)
    rdf.append(x)
    rdf=rdf+list(compose_array)
    rdf.append(sum(compose_array))

  n_ele=len(compose_array)
  rdf=np.array(rdf)  
  rdf=rdf.reshape(-1,n_ele+2)
  return rdf  #x ellist total


def count_neishell_compose_for_givenid(atoms,atm_id,dis_range):

  nei_ids_io=get_neishell_indices(atoms,atm_id,dis_range)
  compose_array=analysis_compose(atoms,nei_ids_io)

  return compose_array

def analysis_compose(atoms,indexs):     
  #indexs=int(indexs)
  eles_tot=np.unique(atoms.get_chemical_symbols())
  #print(eles_tot)
  if len(indexs)>0:
    eles_inlist=atoms[indexs].get_chemical_symbols()
    eles_inlist=np.array(eles_inlist)
    n_ele_list=[]
    for ele in eles_tot:
      ele_idxs=np.where(eles_inlist==ele)[0]
      n_ele=len(ele_idxs)
      n_ele_list.append(n_ele)
    #print(n_ele_list)
  else:
    n_ele_list=list(np.zeros(len(eles_tot)).astype(int)) 
    #print(n_ele_list) 

  n_ele_array=np.array(n_ele_list)
  return n_ele_array

def label_cordpref_info(atoms):
  #work onlyfor fcc
  cutoff=2.85
  dis_range=[1.85,2.85]
  atoms_cordpref={}
  atoms_cordpref['uniq_ele']=list(np.unique(atoms.get_chemical_symbols()))
  atoms_cordpref['uniq_ele_n']=[]
  for ele in atoms_cordpref['uniq_ele']:
    indexs=[atom.index for atom in atoms if atom.symbol == ele]  
    atoms_cordpref['uniq_ele_n'].append(len(indexs))  

  for atm_id in range(len(atoms)):
    atoms_cordpref[atm_id]={}
    compose_array=count_neishell_compose_for_givenid(atoms,atm_id,dis_range)
    CN=sum(compose_array)  
    if CN==6:
      pos_info=['surf','c']
    elif CN==7:
      pos_info=['surf','e']
    elif CN==8:
      pos_info=['surf','111t']
    elif CN==9:
      pos_info=['surf','100t']
    elif CN==12:
      pos_info=['bulk']
    else:
      pos_info=['unrecognised']
    
    atoms_cordpref[atm_id]['ele']=atoms.get_chemical_symbols()[atm_id]
    atoms_cordpref[atm_id]['pos_info']=pos_info
    atoms_cordpref[atm_id]['CN']=CN
    atoms_cordpref[atm_id]['cord_by_ele']=list(compose_array)
    
    print("\rcomplete percentage:{0}% ".format(atm_id*100/len(atoms)), end="", flush=True)
  return atoms_cordpref

def label_cordpref_info_for_Ar(atoms_Vae, good_dGmax_thresh):
  #work onlyfor fcc
  atoms=atoms_Vae
  #cutoff=3.1
  cutoff=2
  #dis_range=[1.85,2.85]
  atoms_cordpref={}
  atoms_cordpref['uniq_ele']=list(np.unique(atoms.get_chemical_symbols()))
  atoms_cordpref['uniq_ele_n']=[]
  for ele in atoms_cordpref['uniq_ele']:
    indexs=[atom.index for atom in atoms if atom.symbol == ele]
    atoms_cordpref['uniq_ele_n'].append(len(indexs))
  indexs_Ar=[atom.index for atom in atoms if atom.symbol == 'Ar']
  indexs_active_Ar=[]
  print("toal Ar sites number are: "+str(len(indexs_Ar)))
  for atm_id in indexs_Ar:
    if atoms.get_initial_charges()[atm_id]<=good_dGmax_thresh:
      indexs_active_Ar.append(atm_id)
  print("Among them, the active (dGmax <"+str(good_dGmax_thresh)+" ) Ar sites number are: "+str(len(indexs_active_Ar)))

  count=1
  atoms_cordpref['id_content']=indexs_active_Ar
  for atm_id in indexs_active_Ar:
    atoms_cordpref[atm_id]={}
    nei_ids_io=get_nei_indices(atoms,atm_id,cutoff)[2:]
    Ar_index_in_nei_ids_io=atoms_cordpref['uniq_ele'].index('Ar')
    compose_array=analysis_compose(atoms,nei_ids_io)
    
    #pos_info,CN=analysis_posinfor_site(atoms,nei_ids_io)
   
    atoms_cordpref[atm_id]['dGmax']=atoms.get_initial_charges()[atm_id]
    atoms_cordpref[atm_id]['ele']=atoms.get_chemical_symbols()[atm_id]
    #atoms_cordpref[atm_id]['pos_info']=pos_info
    #atoms_cordpref[atm_id]['CN']=CN
    atoms_cordpref[atm_id]['nei_ids_io']=nei_ids_io
    atoms_cordpref[atm_id]['cord_by_ele']=list(compose_array)

    print("\rcomplete percentage:{0}% ".format(count*100/len(indexs_active_Ar)), end="", flush=True)
    count=count+1
  return atoms_cordpref

def analysis_posinfor_mat(atoms,atm_id,char_length,first_Ar_id):
  set_cutoff=char_length
  nei_indices=get_nei_indices(atoms,atm_id,set_cutoff)
  nei_indices=nei_indices[2:]

  #print(nei_indices)
  nei_indices_no_ar=[]
  for idx in nei_indices:  
    if idx < first_Ar_id:
      nei_indices_no_ar.append(idx)
  nei_indices=nei_indices_no_ar
  #print(nei_indices)
  CN=len(nei_indices)
  if CN==6:
    pos_info=['surf_c']
  elif CN==7:
    pos_info=['surf_e']
  elif CN==8:
    pos_info=['surf_111t']
  elif CN==9:
    pos_info=['surf_100t']
  elif CN==12:
    pos_info=['bulk']
  else:
    pos_info=['unrecognised']  
  return pos_info,CN

def analysis_posinfor_site(atoms,nei_ids_io,first_Ar_id):
  atoms_nl=atoms[nei_ids_io]
  n_tot=len(atoms_nl)
  indexs_Ar=[atom.index for atom in atoms_nl if atom.symbol == 'Ar']
  n_Ar=len(indexs_Ar)
  CN=n_tot-n_Ar
  nei_ids_mat=[]
  for atm_id in nei_ids_io:
    if atm_id < first_Ar_id:
      nei_ids_mat.append(atm_id)
  #nei_ids_mat: the neibor id that is metal/cat
  pos_info_coll_nei_mat=[]
  nei_eles=atoms[nei_ids_mat].get_chemical_symbols()
  for atm_id in nei_ids_mat:
      char_length=2.85
      pos_info_nei_mat,CN_mat=analysis_posinfor_mat(atoms,atm_id,char_length,first_Ar_id)
      pos_info_coll_nei_mat.append(pos_info_nei_mat)
      
  pos_info=[]
  if CN==1:
    pos_info.append('top')
    pos_info.append(pos_info_coll_nei_mat[0])
  elif CN==2:
    pos_info.append('bri')
    if 'surf_c' in pos_info_coll_nei_mat:
      pos_info.append('c')  
    elif 'surf_e' in pos_info_coll_nei_mat:
      pos_info.append('e')
    else:
      pos_info.append(pos_info_nei_mat[0])
  elif CN==3:
    pos_info.append('3f')
    if 'surf_c' in pos_info_coll_nei_mat:
      pos_info.append('c')
    elif 'surf_e' in pos_info_coll_nei_mat:
      pos_info.append('e')
    else:
      pos_info.append(pos_info_nei_mat[0])
  else:
    pos_info.append('unknown')
  return pos_info,CN,nei_eles

def cal_surf_dis_pref(atoms_cordpref):
  uniq_ele=atoms_cordpref['uniq_ele']
  uniq_ele_n=atoms_cordpref['uniq_ele_n']
  natom=sum(uniq_ele_n)
  surf_ids=[]
  surf_uniq_ele_n=np.zeros(len(uniq_ele))
  for atm_id in range(natom):
    if atoms_cordpref[atm_id]['pos_info'][0]=='surf':
      surf_ids.append(atm_id)
      ele_id=uniq_ele.index(atoms_cordpref[atm_id]['ele'])
      surf_uniq_ele_n[ele_id]=surf_uniq_ele_n[ele_id]+1
  uniq_ele_n=np.array(uniq_ele_n)
  surf_frac=surf_uniq_ele_n/np.sum(surf_uniq_ele_n)
  surf_pref=np.divide(surf_frac,uniq_ele_n)
  surf_pref=surf_pref/np.sum(surf_pref)

  return surf_frac,surf_pref

def output_cord_ele_frac_by_atm_ids(atoms_cordpref,atm_ids):
  count=1
  atm_ids=np.array(atm_ids).astype(int)
  atm_ids=list(atm_ids)

  #atm_ids=int(atm_ids)
  uniq_ele=atoms_cordpref['uniq_ele']
  for atm_id in atm_ids:
   #----------------------count nei ele number -----------#
    cord_by_ele=atoms_cordpref[atm_id]['cord_by_ele']
    if count==1:
      cord_ele_tot=np.array(cord_by_ele)
    else:
      cord_ele_tot=cord_ele_tot+np.array(cord_by_ele)
    count=count+1
  cord_ele_frac=cord_ele_tot/np.sum(cord_ele_tot)
  uniq_ele_n=np.array(atoms_cordpref['uniq_ele_n'])
  cord_ele_pref=np.divide(cord_ele_frac,uniq_ele_n)
  cord_ele_pref=cord_ele_pref/np.sum(cord_ele_pref)

  return cord_ele_frac,cord_ele_pref


def gen_atoms_for_a_given_atm_id(atoms,atm_id, set_cutoff,whether_print=False):
  print("gen_atoms for atm_id of "+str(atm_id)+", with cutoff of "+str(set_cutoff))
  cutOff = neighborlist.natural_cutoffs(atoms)
  for i in range(len(cutOff)):
    cutOff[i]=set_cutoff/2
  nl = neighborlist.NeighborList(cutOff, skin=0, bothways=True)
  nl.update(atoms)
  indices, offsets = nl.get_neighbors(atm_id)
  #----------------------mod 230111 ----------------------------------
  atoms_T=atoms.copy()
  poss_T=atoms_T.get_positions()
  cell=atoms_T.get_cell()
  for atm_jd in range(len(atoms_T)):
    if atm_jd in indices:
      #print("====for id of"+str(atm_jd))
      offset_idx=list(indices).index(atm_jd)
      offset=offsets[offset_idx]
      #print(offset)
      poss_T[atm_jd]=poss_T[atm_jd]+offset[0]*cell[0]+offset[1]*cell[1]+offset[2]*cell[2]
  atoms_T.set_positions(poss_T)
  #----------------------------------------------------------------
  atm_nl=sort_the_indices(atoms_T,atm_id,indices)
  neibor_id_list_of_atm_id=atm_nl[:,0]
  #print(neibor_id_list_of_atm_id)
  ll=[]
  for i in range(len(neibor_id_list_of_atm_id)):
      ll.append(neibor_id_list_of_atm_id[i].astype('int64'))
  atoms_for_atm_id=atoms_T[ll]
  if whether_print==True:
      a=os.path.exists("atoms_nl_coll")
      if a==False:
        b=os.getcwd()
        os.mkdir(b+"/atoms_nl_coll")
      write('./atoms_nl_coll/atom_nl_id'+str(atm_id+1)+'.cif',atoms_for_atm_id)
  return atoms_for_atm_id




def gen_atoms_nl_4_atm_ids(atoms,atm_ids, set_cutoff,whether_print=False):

  #this funcitonal is created at 230515. to make the espression more clear.
  print("now using gen_atoms_for_atm_ids")

  cutOff = neighborlist.natural_cutoffs(atoms)
  #print(0)
  for i in range(len(cutOff)):
    cutOff[i]=set_cutoff/2
  #print(1)
  nl = neighborlist.NeighborList(cutOff, skin=0, bothways=True)

  #print(2)
  nl.update(atoms)
  atoms_nls={}

  #print(3)
  #print("the progess for the gen_atoms_for_each_atom: ", end=" ")


  if atm_ids=='all':
    atm_ids=list(range(len(atoms)))

  count=0
  for atm_id in atm_ids:
    count=count+1

    print("\rcomplete percentage:{0}% ".format(count*100/len(atm_ids)), end="", flush=True)

    atoms_nl=make_atoms_nl(atoms,nl,atm_id)
    atoms_nls[atm_id]=atoms_nl
    

    if whether_print==True:
      a=os.path.exists("atoms_nl_coll")
      if a==False:
        b=os.getcwd()
        os.mkdir(b+"/atoms_nl_coll")
      write('./atoms_nl_coll/atom_nl_id'+str(atm_id+1)+'.cif',atoms_nls[atm_id])


  return atoms_nls,nl



def make_atoms_nl (atoms,nl,atm_id):

  indices, offsets = nl.get_neighbors(atm_id)
  the_id=np.where(indices == atm_id)[0]
  #plus_offset=offsets[plus_id]  
  #indices=np.insert(indices, 0, np.array(atm_id) , axis=0)
  #offsets=np.insert(offsets, 0, plus_offset , axis=0)


  cell=atoms.get_cell()
 
  #print(indices)  
 
  p1=atoms[indices].get_positions()
  D,Dlen=get_distances(p1, cell=cell, pbc='True')

  
  sorted_idx=np.argsort(Dlen[the_id,:])[0]
 
  #print("sorted_idx") 
  #print(sorted_idx)

  indices = indices[sorted_idx]
  offsets = offsets[sorted_idx]

  atoms_nl=atoms[indices]

  cors=atoms_nl.get_positions()

  for i in range(len(atoms_nl)):
    
    offset=offsets[i]
    cors[i]=cors[i]+offset[0]*cell[0]+offset[1]*cell[1]+offset[2]*cell[2]

  atoms_nl.set_positions(cors)


  return atoms_nl



def gen_atoms_for_each_atom(atoms,set_cutoff,whether_print=False):
  #this funcitonal will generate lists of the neighbor atoms (within radii of set_cutoff) for all the atoms in file atoms. And the list will be sorted.
  print("now using gen_atoms_for_each_atom")
  cutOff = neighborlist.natural_cutoffs(atoms)
  print(0)
  for i in range(len(cutOff)):
    cutOff[i]=set_cutoff/2   
  print(1)
  nl = neighborlist.NeighborList(cutOff, skin=0, bothways=True)
  print(2)
  nl.update(atoms)
  atoms_for_each_atom={}
  atm_nl={}
  print(3)
  print("the progess for the gen_atoms_for_each_atom: ", end=" ")
  for atm_id in range(len(atoms)):
    indices, offsets = nl.get_neighbors(atm_id)
    atm_nl[atm_id]=sort_the_indices(atoms,atm_id,indices)
    neibor_id_list_of_atm_id=atm_nl[atm_id][:,0]
    ll=[]
    print("\rcomplete percentage:{0}% ".format(atm_id*100/len(atoms)), end="", flush=True)
    for i in range(len(neibor_id_list_of_atm_id)):
      ll.append(neibor_id_list_of_atm_id[i].astype('int64'))
    atoms_for_each_atom[atm_id]=atoms[ll]
    if whether_print==True:
      a=os.path.exists("atoms_nl_coll")
      if a==False:
        b=os.getcwd()
        os.mkdir(b+"/atoms_nl_coll")
      write('./atoms_nl_coll/atom_nl_id'+str(atm_id+1)+'.cif',atoms_for_each_atom[atm_id])
  return atoms_for_each_atom




def inside_the_box(cor,box):
  #box=[veca,vecb,vecc,corO]
  corO=box[3]
  relative_cor=cor-corO
  #move to origin
  cell=box[0:3]
  cell=np.mat(cell)
  scaled_cor=relative_cor*np.linalg.inv(cell)
  #print(scaled_cor)

  if np.max(scaled_cor)<0.999 and np.min(scaled_cor)>=0:
    result=1
  else:
    result=0
  return result



def gen_atoms_for_given_box(atoms,box,whether_print=False):
  corO=box[3]
  cell=box[0:3]
  atm_ids_in_box=[]
  for atm_jd in range(len(atoms)):
    if inside_the_box(atoms.get_positions()[atm_jd],box):
      atm_ids_in_box.append(atm_jd)
  atoms_in_box=atoms[atm_ids_in_box]
  atoms_in_box.translate(-corO)
  atoms_in_box.set_cell(cell)
  if whether_print==True:
    write('atom_nl_box.cif',atoms_in_box)
  return atoms_in_box


def sort_the_indices(atoms,i,indices):
  A=[]
  for index in indices:
    dis=atoms.get_distance(i,index,mic=True)
    A=np.append(A,[index,dis])
  A=A.reshape(-1,2)
  A=A[np.argsort(A[:,1])]
  for j in range(len(A)):
    A[j,0]=int(A[j,0])
  neibor_indices_sort_by_dis=A
  return neibor_indices_sort_by_dis

def gen_e_vec(atoms,mode=0):
  #------------------------beginiing of generation of x_vec and nei_cl----------------------
  mode_dependent_list=[[2,3],[2,4],[3,4],[2,5],[3,5],[4,5],[3,2],[4,2],[4,3],[5,2],[5,3],[5,4]]
  #mode_dependent_list=[[3,2],[4,2],[4,3],[5,2],[5,3],[5,4]]

  mode_dependent_list=np.array(mode_dependent_list)
  e_vec=[]
  if len(atoms)==1 or len(atoms)==2:
    ei1=[1,0,0]
    ei2=[0,1,0]
    ei3=[0,0,1]
  elif len(atoms)==3:
    Ria=atoms.get_distance(0,2,mic=True,vector=True)   #notice neibor_id_list_of_atm_id[0] is itself !
    ei1=Ria/np.linalg.norm(Ria)
    Rib=[0,1,1]
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  else:
    if mode<=6:
      id1=min([mode_dependent_list[mode,0],len(atoms)-2])
      id2=min([mode_dependent_list[mode,1],len(atoms)-1])
    elif mode>6 and mode<=12:
      id1=min([mode_dependent_list[mode,0],len(atoms)-1])
      id2=min([mode_dependent_list[mode,1],len(atoms)-2])
    Ria=atoms.get_distance(0,id1,mic=True,vector=True)
    ei1=Ria/np.linalg.norm(Ria)
    Rib=atoms.get_distance(0,id2,mic=True,vector=True)
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)

  e_vec=[ei1,ei2,ei3]
  return e_vec

def gen_e_vec_beta(atoms):
  #------------------------beginiing of generation of x_vec and nei_cl----------------------
  e_vec=[]
  if len(atoms)==1 or len(atoms)==2:
    ei1=[1,0,0]
    ei2=[0,1,0]
    ei3=[0,0,1]
  elif len(atoms)==3:
    Ria=atoms.get_distance(0,2,mic=True,vector=True)   #notice neibor_id_list_of_atm_id[0] is itself !
    ei1=Ria/np.linalg.norm(Ria)
    Rib=[0,1,1]
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  else:
    Ria=atoms.get_distance(0,3,mic=True,vector=True)
    ei1=Ria/np.linalg.norm(Ria)
    Rib=atoms.get_distance(0,2,mic=True,vector=True)
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  e_vec=[ei1,ei2,ei3]
  return e_vec


def gen_Dij(atoms,mode=0):
  e_vec=gen_e_vec(atoms,mode)
  Dij=[]
  symbols=[]
  for atm_id in range(len(atoms)):
    if atm_id<=1:
      Dij=np.append(Dij,[1,0,0,0])
      symbols.append(atoms.get_chemical_symbols()[atm_id])
    else:
      #Rij=atoms.get_distance(0,atm_id,mic=True)
      Rij_vec=atoms.get_distance(0,atm_id,mic=True,vector=True)


      Rij_vec=np.mat(Rij_vec)
      #how to measure cor?
      #[ei1,ei2,ei3]*[xij;yij;zij]=[Rij_vec]'  3x3  3x1  3x1
      #=>[xij;yij;zij]=[x_vec,y_vec,z_vec]^-1*[Rij_vec]'
      xij_yij_zij=Rij_vec*np.linalg.inv(e_vec)
      xij_yij_zij=xij_yij_zij.reshape(-1,3) 
      
      xij= xij_yij_zij[0,0]
      yij= xij_yij_zij[0,1]
      zij= xij_yij_zij[0,2]
      Dij=np.append(Dij,[1,xij,yij,zij])
      symbols.append(atoms.get_chemical_symbols()[atm_id])
  #if len(Dij)>=1:
  Dij=Dij.reshape(-1,4)
  return Dij,symbols,e_vec

def print_atoms_under_all_modes(atoms,prefix):
  for mode in range(6):
    e_vec=gen_e_vec(atoms,mode) 
    rotted_atoms=atoms_rotted_by_its_e_evc(atoms,e_vec)
    write(prefix+"_mode"+str(mode)+'.xyz',rotted_atoms)
  return

def score_the_similarity(atomsa,atomsb):
  mode_scores=[]
  for mode1 in range(12):
    for mode2 in range(12):
      Dija,symbolsa,e_vec_a=gen_Dij(atomsa,mode1)
      Dijb,symbolsb,e_vec_b=gen_Dij(atomsb,mode2)
      score=score_the_similarity_detailed(Dija,symbolsa,Dijb,symbolsb)
      mode_scores=np.append(mode_scores,[mode1,mode2,score])
  mode_scores=mode_scores.reshape(-1,3)
  mode_scores=mode_scores[np.argsort(mode_scores[:,2])]
  score=mode_scores[0,2]
  #print("mode_scores")
  #print(mode_scores)
  mode_of_a=mode_scores[0,0].astype('int64')
  mode_of_b=mode_scores[0,1].astype('int64')

  return score,mode_of_a,mode_of_b


def score_the_similarity_detailed(Dija,symbolsa,Dijb,symbolsb):

  if len(Dija)<len(Dijb) or symbolsa[0]!=symbolsb[0]:
    #print("length of Dij is not match or central element dont match , skip")
    score=100
  else:
    AA=[]
    for atmb_id in range(len(Dijb)):
      specieb=symbolsb[atmb_id]
      A=[]
      for atma_id in range(len(Dija)): 
        speciea=symbolsa[atma_id]
        if speciea==specieb:
          dis=np.linalg.norm(Dija[atma_id]-Dijb[atmb_id])
          dis=dis**2
          A=np.append(A,[atma_id,dis])
      if len(A)==0:
        smallest_dis_of_atma_ids_to_atmb_id=100
      else:
        A=A.reshape(-1,2)
        A=A[np.argsort(A[:,1])]
        smallest_dis_of_atma_ids_to_atmb_id=A[0,1]
      AA=np.append(AA,[atmb_id,smallest_dis_of_atma_ids_to_atmb_id])
    AA=AA.reshape(-1,2)
    AA=AA[np.argsort(AA[:,1])]
    score=AA[len(AA)-1,1]


  return score




def merge_based_on_individual_similar_part_no_del_atomsb_nl(atomsa,atomsa_nl,atomsb,atomsb_nl,mode_of_a=0,mode_of_b=0):
  #score,mode_of_a,mode_of_b=score_the_similarity(atomsa_nl,atomsb_nl)
  print("merging")
  print("mode_of_a,mode_of_b")
  print(mode_of_a,mode_of_b)
  e_vec_a=gen_e_vec(atomsa_nl,mode_of_a)
  e_vec_b=gen_e_vec(atomsb_nl,mode_of_b)

  impointb=atomsb_nl.get_positions()[0]
  impoint_ref=atomsa_nl.get_positions()[0]
  e_vec_ref=e_vec_a
  rotrance_atomsb=rotrance_atoms_to_ref(atomsb,impointb,impoint_ref,e_vec_b,e_vec_ref)
  new_atoms=merge_atoms(atomsa,rotrance_atomsb)

  return new_atoms

def merge_based_on_individual_similar_part(atomsa,atomsa_nl,atomsb,atomsb_nl,mode_of_a=0,mode_of_b=0,del_cutoff=0.2):
  #score,mode_of_a,mode_of_b=score_the_similarity(atomsa_nl,atomsb_nl)
  print("merging")
  print("mode_of_a,mode_of_b")
  print(mode_of_a,mode_of_b)
  e_vec_a=gen_e_vec(atomsa_nl,mode_of_a)
  e_vec_b=gen_e_vec(atomsb_nl,mode_of_b)

  impointb=atomsb_nl.get_positions()[0]
  impoint_ref=atomsa_nl.get_positions()[0]
  e_vec_ref=e_vec_a
  rotrance_atomsb=rotrance_atoms_to_ref(atomsb,impointb,impoint_ref,e_vec_b,e_vec_ref)
  new_atoms=merge_atoms(atomsa,rotrance_atomsb)
  get_duplicate_atoms(new_atoms, del_cutoff, delete=True)

  return new_atoms


#def substract_atoms(atomsa,atomsb):
#  for atmb_id in range(len(atomsb)):
#    A=[]
#    for atma_id in range(len(atomsa)):
#      dis=np.linalg.norm(atomsa.get_positions()[atma_id]-atomsb.get_positions()[atmb_id])
#      A=np.append(A,[atma_id,dis])
#    A=A.reshape(-1,2)
#    A=A[np.argsort(A[:,1])]
#    del_atm_id=A[0,0]      


  return new_atoms

def merge_atoms(atomsa,atomsb):
  symbolesa=atomsa.get_chemical_symbols()
  symbolesb=atomsb.get_chemical_symbols()
  symboles=symbolesa+symbolesb

  positions=np.append(atomsa.get_positions(),atomsb.get_positions())
  positions=positions.reshape(-1,3)
  new_atoms=Atoms(symboles,positions,cell=atomsa.get_cell(),pbc=True)
  return new_atoms

def merge_atoms_and_del_dup(atomsa,atomsb,del_cutoff=0.1):
  symbolesa=atomsa.get_chemical_symbols()
  symbolesb=atomsb.get_chemical_symbols()
  symboles=symbolesa+symbolesb

  positions=np.append(atomsa.get_positions(),atomsb.get_positions())
  positions=positions.reshape(-1,3)
  new_atoms=Atoms(symboles,positions,cell=atomsa.get_cell(),pbc=True)
  get_duplicate_atoms(new_atoms, del_cutoff, delete=True)

  return new_atoms

def opt_too_close_atoms_to_calculatable_one(atoms,dtoler=0.3,ds=0.2):
  
  new_atoms=atoms.copy()
  d=get_the_min_dis_in_atoms(atoms)
  count=0
  print("iter "+str(count)+', we got the min dis to be '+str(d))
  while d < dtoler:
    print("d is larger than the dtoler of "+str(dtoler))
    i_j_dis=[]
    for atm_id in range(len(atoms)):
      for atm_jd in range(atm_id+1,len(atoms)):
        dis=atoms.get_distance(atm_id,atm_jd,mic=True)
        i_j_dis=np.append(i_j_dis,[atm_id,atm_jd,dis]) 
    i_j_dis=i_j_dis.reshape(-1,3)
    arg_idx=np.argsort(i_j_dis[:,2])
    i_j_dis=i_j_dis[arg_idx]
    i_j_dis_chose=i_j_dis[0]
    i=int(i_j_dis_chose[0])
    j=int(i_j_dis_chose[1])
    print("expand the dis btewenn "+str(i)+' and '+str(j))
    vec_ij=atoms.get_distance(i, j, mic=True, vector=True)
    vec_ij=vec_ij/np.linalg.norm(vec_ij)
    positions=atoms.get_positions()
    positions[i]=positions[i]-vec_ij*ds
    positions[j]=positions[j]+vec_ij*ds
    atoms.set_positions(positions)
    d=get_the_min_dis_in_atoms(atoms)
    count=count+1
    print("iter "+str(count)+', we got the min dis to be '+str(d))
  print("iter reches dtoler at iter "+str(count))
  new_atoms=atoms
  #d=get_the_min_dis_in_atoms(new_atoms)
  #print("final dis mmin is "+str(d))
  return new_atoms

def rotrance_atoms_to_ref(atoms,impoint,impoint_ref,e_vec,e_vec_ref):
  e_vec=np.mat(e_vec)
  e_vec_ref=np.mat(e_vec_ref)
  e_vec=e_vec.reshape(-1,3)
  e_vec_ref=e_vec_ref.reshape(-1,3)
  #e_vec=np.mat(e_vec)
  #e_vec_ref=np.mat(e_vec_ref)
  #base on:
  #(atoms.get_positions()-impoint)*np.linalg.inv(e_vec)*e_vac_ref+impoint_ref
  rotrance_cor=(atoms.get_positions()-impoint)*np.linalg.inv(e_vec)*e_vec_ref+impoint_ref
  symboles=atoms.get_chemical_symbols()
  rotrance_atoms=Atoms(symboles,rotrance_cor,cell=atoms.get_cell(),pbc=True)

  return rotrance_atoms



def cal_the_tranM(e_vec_a,e_vec_b):
  e_vec_a=e_vec_a.reshape(-1,3)
  e_vec_b=e_vec_b.reshape(-1,3)
  e_vec_a=np.mat(e_vec_a)
  e_vec_b=np.mat(e_vec_b)

  e_vec_a_inv=np.linalg.inv(e_vec_a)
  M_trans=e_vec_a_inv*e_vec_b

  return M_trans

def switch_symbol_to_HHeli(atoms):
  symbols=atoms.get_chemical_symbols()
  symbols_switch=[]
  for i in range(len(atoms)):
    #if symbols[i]=='H':
    #  symbols_switch.append('Li')
    #if symbols[i]=='O':
    #  symbols_switch.append('Li')
    #if symbols[i]=='C':
    #  symbols_switch.append('Li')
    #if symbols[i]=='Co':
    #  symbols_switch.append('Li')
    #if symbols[i]=='N':
    #  symbols_switch.append('Li')
    #if symbols[i]=='Bi':
    #  symbols_switch.append('Li')
    symbols_switch.append('Nb')
  #print(symbols_switch)
  positions=atoms.get_positions()
  new_atoms=Atoms(symbols,positions,cell=atoms.get_cell(),pbc=True)
  new_atoms.set_chemical_symbols(symbols_switch)
  return new_atoms

def get_the_max_dis_in_atoms(atoms,mic_value=True):
  A=[]
  for atm_id in range(len(atoms)):
    for atm_jd in range(atm_id+1,len(atoms)):
      A.append(atoms.get_distance(atm_id,atm_jd,mic=mic_value))
  max_dis=max(A)
  return max_dis

def get_the_min_dis_in_atoms(atoms,mic_value=True):
  A=[]
  for atm_id in range(len(atoms)):
    for atm_jd in range(atm_id+1,len(atoms)):
      A.append(atoms.get_distance(atm_id,atm_jd,mic=mic_value))
  min_dis=min(A)
  return min_dis

def get_the_min_dis_in_atoms_for_atm_id(atoms,atm_id):
  A=[]
  for atm_jd in range(len(atoms)):
    if atm_jd!=atm_id:
      A.append(atoms.get_distance(atm_id,atm_jd,mic=True))
  min_dis=min(A)
  return min_dis

def get_the_max_dis_in_atoms_for_atm_id(atoms,atm_id):
  A=[]
  for atm_jd in range(len(atoms)):
    if atm_jd!=atm_id:
      A.append(atoms.get_distance(atm_id,atm_jd,mic=True))
  max_dis=min(A)
  return max_dis

def largest_distance_in_a_given_dir(cor,dir_vec):
  #cor should be in fmt of np.array
  A=[]
  atm_id=0
  for cor_i in cor:
    dis=np.dot(cor_i,dir_vec) 
    A=np.append(A,[atm_id,dis])
    atm_id=atm_id+1
  A=A.reshape(-1,2)
  A=A[np.argsort(A[:,1])]
  atm_id_of_max_dis=A[-1,0]
  max_dis=A[-1,1]
  atm_id_of_max_dis_neg=A[0,0]
  max_dis_neg=A[0,1]
  return max_dis,max_dis_neg


def opt_cell_of_xy_right_angle_version(atoms):
  #used in non-pbc, but atoms should already get cell
  edge_length=1.5
  center_cor=np.mean(atoms.get_positions(),axis=0)
  cor=atoms.get_positions()
  cor=cor-center_cor
  i_S=[]
  for i in range(0,18):
    thita=np.pi/36*i    # 0 to pi/2
    x1=np.cos(thita)
    y1=np.sin(thita)
    x2=np.cos(thita+np.pi/2)
    y2=np.sin(thita+np.pi/2)
    vec1=[x1,y1,0]
    vec1=np.array(vec1)
    vec2=[x2,y2,0]
    vec2=np.array(vec2)
    max_dis1,max_dis_neg1=largest_distance_in_a_given_dir(cor,vec1)
    max_dis2,max_dis_neg2=largest_distance_in_a_given_dir(cor,vec2)
    a=max_dis1-max_dis_neg1
    b=max_dis2-max_dis_neg2
    S=a*b
    i_S=np.append(i_S,[i,S])
  i_S=i_S.reshape(-1,2)
  i_S=i_S[np.argsort(i_S[:,1])] 
  i_of_Smin=i_S[0,0]
  i=i_of_Smin
  #-----------------------------------------------------------
  thita=np.pi/36*i
  x1=np.cos(thita)
  y1=np.sin(thita)
  x2=np.cos(thita+np.pi/2)
  y2=np.sin(thita+np.pi/2)
  vec1=[x1,y1,0]
  vec1=np.array(vec1)
  vec2=[x2,y2,0]
  vec2=np.array(vec2)
  vec3=[0,0,1]
  vec3=np.array(vec3)
  max_dis1,max_dis_neg1=largest_distance_in_a_given_dir(cor,vec1)
  max_dis2,max_dis_neg2=largest_distance_in_a_given_dir(cor,vec2)
  max_dis3,max_dis_neg3=largest_distance_in_a_given_dir(cor,vec3)
  a=max_dis1-max_dis_neg1
  b=max_dis2-max_dis_neg2
  c=max_dis3-max_dis_neg3 
  #---------------------------------------------------------------
  cor=cor+(a-max_dis1+edge_length)*vec1+(b-max_dis2+edge_length)*vec2+(edge_length-max_dis_neg3)*vec3  #put it into box center
  atoms.set_positions(cor)
  cella=(a+2*edge_length)*vec1
  cellb=(b+2*edge_length)*vec2
  cellc=[0,0,c+20]
  cellc=np.array(cellc)
  cell=np.vstack([cella,cellb,cellc]) 
  atoms.set_cell(cell) 
    
  return atoms  

def opt_vac_thickness(atoms):
  print("=====================opt_vac_thickness()=================")
  edge_length=1.5
  c=atoms.get_cell()[2,2]
  ref_thickness=6  #(thickness of vac)
  print("opt_vac_thickness(): try to adjust the thickness of vac to",ref_thickness)
  #print("opt_vac_thickness(): orignal lenght of c in cell is "+str(c))
  corz_coll=[]
  corz_coll=atoms.get_positions()[:,2]
  corz_coll=np.sort(corz_coll)
  corz_coll=list(corz_coll)
  corz_coll.append(corz_coll[0]+c)
  #print("corz_coll")
  #print(corz_coll)
  z_diszdown=[]  #the z and its dis inz for the atoms downwards
  for i in range(len(corz_coll)):
    diszdown=corz_coll[i]-corz_coll[i-1]
    z_diszdown=np.append(z_diszdown,[corz_coll[i],diszdown])
  z_diszdown=z_diszdown.reshape(-1,2)
  z_diszdown=z_diszdown[np.argsort(z_diszdown[:,1])]
  ref_vec_thick=z_diszdown[-1,1]
  print("opt_vac_thickness(): first, find the referece thickness (largest distance across z-pbc) is "+str(ref_vec_thick),end="")
  ref_z=z_diszdown[-1,0] # above which should minus ref_vec_thick-9
  print(", with the location of vac is at "+str(ref_z)+" and "+str(ref_z+ref_vec_thick))
  dzvalue_of_all=c-ref_z+1.5  #up shift this value of all to keep all vac layer along z=0

  if ref_vec_thick >= ref_thickness:
    print("opt_vac_thickness(): the referece thickness is large than expected thickness ("+str(ref_thickness)+") , should adjust iti. ")
    positions=atoms.get_positions()
    shift_in_z=ref_thickness-ref_vec_thick
    print("opt_vac_thickness(): we adjust it by down shifting of "+str(shift_in_z)+"(the difference between referece and expected thickness)",end="")
    print(", for all the atoms has z cor larger than "+str(ref_z-1))
    for atm_id in range(len(atoms)):
      if positions[atm_id,2]>=ref_z-1:
        print("opt_vac_thickness(): from "+str(positions[atm_id,2]),end="")
        positions[atm_id,2]=positions[atm_id,2]+shift_in_z
        print("to "+str(positions[atm_id,2]))
    atoms.set_positions(positions)
    cell=atoms.get_cell()
    print("opt_vac_thickness(): finally, shorten OC, from "+str(cell[2,2]),end="")
    cell[2,2]=cell[2,2]+shift_in_z
    print("to "+str(cell[2,2]))
    atoms.set_cell(cell)
  else:
    print("opt_vac_thickness(): referece thickness is smaller than expected thickness, no need to adjust")
  #-------------------------------------------put the vac at the top of the box
  print("opt_vac_thickness(): now put the vac at the top of the box")
  positions=atoms.get_positions()
  c=atoms.get_cell()[2,2]
  print("opt_vac_thickness(): the new c is"+str(c))
  print("opt_vac_thickness(): up shift all the atoms by "+str(dzvalue_of_all))
  for atm_id in range(len(atoms)):
    positions[atm_id,2]=positions[atm_id,2]+dzvalue_of_all
    if positions[atm_id,2]>c:
      print("z value of "+str(atm_id)+"is "+str(positions[atm_id,2])+", it is larger than c, downshift to the bottom box")
      positions[atm_id,2]=positions[atm_id,2]-c
  atoms.set_positions(positions)

  print("=====================opt_vac_thickness().end=================")
  
  return atoms


def smallest_cell_make(atoms):
  edge_length=1.5
  #thie script try to make an appropriate cell based on the vilumn of the moleclues
  center_cor=np.mean(atoms.get_positions(),axis=0)
  center_corid=get_center_atm_id(atoms)
  #print("center_corid")
  #print(center_corid)
  #1.vec a make
  #write('atom_Test_before.xyz',atoms)
  #print("1. make vec a based on the larest distance of the the cell")
  A=[]
  for atm_id in range(len(atoms)):
    dis=np.linalg.norm(atoms.get_positions()[atm_id]-center_cor)
    A=np.append(A,[atm_id,dis])
  A=A.reshape(-1,2)
  A=A[np.argsort(A[:,1])]
  farest_atm_id=A[-1,0]
  #print("farest_atm_id for x")
  #print(farest_atm_id)
  farest_atm_id=int(farest_atm_id)
  farest_dis=A[-1,1]
  len_a=farest_dis+edge_length
  e_veca=(atoms.get_positions()[farest_atm_id]-center_cor)/farest_dis
  #2.vec b make
  print("2. make vec a based on the larest distance to axis vec a")
  A=[]
  for atm_id in range(len(atoms)):
    vec_b_test=atoms.get_positions()[atm_id]-center_cor
    dis=np.linalg.norm(np.cross(vec_b_test,e_veca))
    A=np.append(A,[atm_id,dis])
  A=A.reshape(-1,2)
  A=A[np.argsort(A[:,1])]
  farest_atm_id=A[-1,0]
  farest_atm_id=int(farest_atm_id)
  farest_dis=A[-1,1]
  len_b=farest_dis+edge_length
  b=atoms.get_positions()[farest_atm_id]-center_cor
  e_vecb=b-np.dot(b,e_veca)*e_veca
  e_vecb=e_vecb/np.linalg.norm(e_vecb)
  #e_vecb=(atoms.get_positions()[farest_atm_id]-center_cor)/farest_dis
  #3. vec c make
  print("3. make vec c. c is a little special. we first get the new axis based on vec_b vec_a")
  e_vecc=np.cross(e_veca,e_vecb)
  M_e_vec=[e_veca,e_vecb,e_vecc]
  M_e_vec=np.mat(M_e_vec)
  M_e_vec=M_e_vec.reshape(-1,3)
  print("then we project the atoms to get the new atoms")
  
  new_positions=(atoms.get_positions()-center_cor)*np.linalg.inv(M_e_vec)
  new_positions=np.array(new_positions)
  #print("new_positions")
  #print(new_positions)
  atoms.set_positions(new_positions)
  #write('atom_Test.xyz',atoms)
  print("M_e_vec is ")
  print(M_e_vec)
  A=[]
  for atm_id in range(len(atoms)):
    corz=atoms.get_positions()[atm_id][2]
    dis=abs(corz)
    A=np.append(A,[atm_id,dis])
  A=A.reshape(-1,2)
  A=A[np.argsort(A[:,1])]
  #print("A")
  #print(A)
  farest_atm_id=A[-1,0]
  farest_atm_id=int(farest_atm_id)
  farest_dis=A[-1,1]
  len_c=farest_dis+edge_length
  #e_vecc=(atoms.get_positions()[farest_atm_id]-center_cor)/farest_dis
  #4.collect e_vecx len_x to get cell
  len_a_part2=abs(min(atoms.get_positions()[:,0]))+edge_length
  len_b_part2=abs(min(atoms.get_positions()[:,1]))+edge_length
  len_c_part2=abs(min(atoms.get_positions()[:,2]))+edge_length
  OA=[len_a+len_a_part2,0,0]
  OB=[0,len_b+len_b_part2,0]
  OC=[0,0,len_c+len_c_part2]
  cell=[OA,OB,OC]
  cell=np.array(cell)
  cell=cell.reshape(-1,3)
  atoms.set_cell(cell)
  #write('atom_Test_final.xyz',atoms)
  print("finally, we have cell of")
  print(cell)
  print("now we can shift the corcenter to cell center")
  
  cell_center=[len_a_part2,len_b_part2,len_c_part2]
  cell_center=np.array(cell_center)
  new_positions=atoms.get_positions()+cell_center
  atoms.set_positions(new_positions)

  return atoms


def get_center_atm_id(atoms):
  center_cor=np.mean(atoms.get_positions(),axis=0)
  A=[]
  for atm_id in range(len(atoms)):
   dis=np.linalg.norm(atoms.get_positions()[atm_id]-center_cor)
   A=np.append(A,[atm_id,dis])
  A=A.reshape(-1,2)
  A=A[np.argsort(A[:,1])]
  center_atm_id=A[0,0]
  return center_atm_id


#def cal_ave_corNs(atoms,nl)
#


#def cut_off_test_for_surf_identifier(atoms):
  
#  set_cutoff= 7

#  cutOff = neighborlist.natural_cutoffs(atoms)
#  for i in range(len(cutOff)):
#    cutOff[i]=set_cutoff/2
#  nl = neighborlist.NeighborList(cutOff, skin=0, bothways=True)
#  nl.update(atoms)

#  atm_ids=random.sample(len(atoms),10)  

#  for 




def identify_the_surface_site(atoms,atm_ids='all'):

  set_cutoff= 9

  nl = ini_nl(atoms,set_cutoff)

  if atm_ids=='all':
    atm_ids=range(len(atoms))


  surf_atm_ids=[]

  count=0

  for atm_id in atm_ids:

    count=count+1

    char_vec=0
    nei_ids=get_nei_indices_v2(atoms,atm_id,set_cutoff,nl)
    nei_ids=np.array(nei_ids)    
    nei_ids=nei_ids[nei_ids != atm_id]
    vecs=[]
    for atm_jd in nei_ids:
      R_vec=atoms.get_distance(atm_id,atm_jd,mic=True,vector=True)
      char_vec=char_vec+R_vec/(np.linalg.norm(R_vec)**1.5)
      vecs=np.append(vecs,R_vec)    

    vecs = vecs.reshape(-1,3)
    most_empty_dir=np.array(-char_vec)
    most_empty_dir_repeat = np.tile(most_empty_dir, (len(vecs), 1)) 
    
   
    #print("vecs= ",vecs)
    #print("most_emtry_dir,",most_empty_dir)  
    
    angles=get_angles(vecs, most_empty_dir_repeat)
    #print(angles)

    min_angle=min(angles)
    #print("min_angle=",min_angle)

    if min_angle < 35:
      if_surf=False
    else:
      if_surf=True

    #print(atm_id," is ",if_surf)
    if if_surf:
      surf_atm_ids.append(atm_id)

    print("\ridentify_the_surface_site complete percentage:{0}% ".format(count*100/len(atm_ids)), end="", flush=True)

  return surf_atm_ids






def check_whether_surface_site_normal_version(atoms,atm_id):
  char_vec=0
  
  nl = ini_nl(atoms,9)

  nei_ids=get_nei_indices_v2(atoms,atm_id,9,nl)
  for atm_jd in nei_ids:
    R_vec=atoms.get_distance(atm_id,atm_jd,mic=True,vector=True)
    char_vec=char_vec+R_vec
  char_vec=char_vec/(len(atoms)-2)
  result="true"
  for atm_jd in nei_ids:
    test_vec=atoms.get_positions()[atm_id]-atoms.get_positions()[0]
    cos_value=np.dot(char_vec,test_vec)/np.linalg.norm(char_vec)/np.linalg.norm(test_vec)
    if cos_value<-0.5:
      result="false"
      break
  return result


def check_whether_surface_site(atoms_nl):
  #center_cor=np.mean(atoms.get_positions()) 
  #print("now doing check_whether_surface_site")
  char_vec=0
  for atm_id in range(2,len(atoms_nl)):
    R_vec=atoms_nl.get_distance(0,atm_id,mic=True,vector=True)
    char_vec=char_vec+R_vec
  char_vec=char_vec/(len(atoms_nl)-2)
  #print("generate the char vector for each atoms (0-->atm_id)")  
  #print(char_vec)
  #char_vec=center_cor-atoms.get_positions()[0]
  result="true"
  for atm_id in range(2,len(atoms_nl)):
    test_vec=atoms_nl.get_positions()[atm_id]-atoms_nl.get_positions()[0]
    cos_value=np.dot(char_vec,test_vec)/np.linalg.norm(char_vec)/np.linalg.norm(test_vec)
    if cos_value<-0.5:
      result="false"
      break
  return result


def atoms_rotted_by_its_e_evc(atoms,e_vec):
  relative_positions=[]
  for atm_id in range(len(atoms)):
    Rij_atm_id=atoms.get_distance(0,atm_id,mic=True,vector=True)
    relative_positions=np.append(relative_positions,Rij_atm_id)
  relative_positions=relative_positions.reshape(-1,3)
  e_vec=np.mat(e_vec)
  relative_positions=np.mat(relative_positions)
  rotted_positions=relative_positions*np.linalg.inv(e_vec)+atoms.get_positions()[0]
  rotted_positions=np.array(rotted_positions)
  #print(rotted_positions)
  rotted_atoms=atoms

  rotted_atoms.set_positions(rotted_positions)

  return rotted_atoms

#-------------------------------------------------------------------------
def rot_to_z_axis(atoms, c_atm_id, side_length):
  atoms=shift_cell_by_center_atom(atoms, c_atm_id, side_length)
  ref_vec=cal_norm_vector_of_center_atom(atoms, c_atm_id)
  new_atoms,M_rot=rotate_atoms_to_let_ref_vec_along_z_vec(atoms,c_atm_id,ref_vec)
  return new_atoms,M_rot

def rot_to_z_axis_surface_version (atoms, c_atm_id, side_length):
  atoms=shift_cell_by_center_atom(atoms, c_atm_id, side_length)
  ref_vec=vect_for_most_sparse_locally_v2(atoms,c_atm_id)
  new_atoms,M_rot=rotate_atoms_to_let_ref_vec_along_z_vec(atoms,c_atm_id,ref_vec)
  return new_atoms,M_rot


def shift_cell_by_center_atom(atoms, c_atm_id, side_length):
  cor=atoms.get_positions()
  print("shift O point to orint point")
  shift_cor=cor-atoms.get_positions()[c_atm_id]+[side_length/2,side_length/2,side_length/2]  
  print("note if in cif, center atom cor will be [0.5 0.5 0.5]")
  atoms.set_positions(shift_cor)
  atoms.set_cell([(side_length, 0, 0), (0, side_length, 0), (0, 0, side_length)])
  return atoms



def cal_norm_vector_of_center_atom(atoms, c_atm_id):
  center_cor=np.mean(atoms.get_positions(),axis=0)
  cor_of_c_atm_id=atoms.get_positions()[c_atm_id]
  print("center corrdinate is "+str(center_cor))
  print("corrdinate of the selected id is: "+str(cor_of_c_atm_id))
  print("their difference gives the ref_vec: ")
  ref_vec=cor_of_c_atm_id-center_cor
  ref_vec=ref_vec/np.linalg.norm(ref_vec)
  print(ref_vec)
  
  #print("now rotate atoms, let ref_vec along z axis")

  return ref_vec

def cal_norm_vector_of_center_atom_surface_version(atoms, c_atm_id):
  #center_cor=np.mean(atoms.get_positions(),axis=0)
  cor_of_c_atm_id=atoms.get_positions()[c_atm_id]
  set_cutoff=3
  atoms_nl=gen_atoms_for_a_given_atm_id(atoms,c_atm_id, set_cutoff)  
  local_center_cor=np.mean(atoms_nl.get_positions()[1:-1],axis=0)
  print("center corrdinate is "+str(local_center_cor))
  print("corrdinate of the selected id is: "+str(cor_of_c_atm_id))
  print("their difference gives the ref_vec: ")
  ref_vec=cor_of_c_atm_id-local_center_cor
  ref_vec=ref_vec/np.linalg.norm(ref_vec)
  print(ref_vec)

  #print("now rotate atoms, let ref_vec along z axis")

  return ref_vec


def rotate_atoms_to_let_ref_vec_along_z_vec(atoms,c_atm_id,ref_vec):
  cor_of_c_atm_id=atoms.get_positions()[c_atm_id]
  print("rotate_atoms_to_let_ref_vec_along_z_vec")
  print("\nstep 1 shift center atom to origin point: ")
  cor=atoms.get_positions()-cor_of_c_atm_id
  print("\nstep 2 : the rotate could base on the rotation Matrix. if elemental vect is shift from v_ea to v_eb, the M_rotation=inv(v_ea)*v_eb and goal_cor=cor*M_rotation")
  print("make a v_ea")
  print("set the vect that is smallerst in ref_vec to be one")
  print("ref_vec:")
  print(ref_vec)
  ref_vecx=[0,0,0]
  #-----------------------------if have time re write this-------------  
  ref_vecx[list(np.abs(ref_vec)).index(min(np.abs(ref_vec)))]=1
  ref_vecx[list(np.abs(ref_vec)).index(max(np.abs(ref_vec)))]=0
  print("\nset smallest index to be 1, and largest to be 0")
  print("smallest is "+ str(list(np.abs(ref_vec)).index(min(np.abs(ref_vec))))+", set that of ref_vecx  to be 1")
  print("largest is "+str(list(np.abs(ref_vec)).index(max(np.abs(ref_vec))))+", set that of ref_vecx to 0")
  S=[0,1,2]
  S1=max(list(np.abs(ref_vec)).index(min(np.abs(ref_vec))),list(np.abs(ref_vec)).index(max(np.abs(ref_vec))))
  S2=min(list(np.abs(ref_vec)).index(min(np.abs(ref_vec))),list(np.abs(ref_vec)).index(max(np.abs(ref_vec))))
  print("orderly pop index of "+ str(S1)+" and "+str(S2))
  S.pop(S1)
  S.pop(S2)
  print("we have S=")
  print(S)
  if ref_vec[S[0]]!=0:
    ref_vecx[S[0]]=-ref_vec[list(np.abs(ref_vec)).index(min(np.abs(ref_vec)))]/ref_vec[S[0]]
  else:
    ref_vecx[S[0]]=0
  #-----------------------end of if hav time---------------------------------------------------
  print("ref_vecx is")
  print(ref_vecx)
  ref_vecx=ref_vecx/np.linalg.norm(ref_vecx)
  print("after norm it becomes")
  print(ref_vecx)
  print("\nget last vec by np.cross")
  ref_vecy=np.cross(ref_vecx,ref_vec)
  ref_vecy=ref_vecy/np.linalg.norm(ref_vecy)
  v_ea=list(ref_vecx)+list(ref_vecy)+list(ref_vec)
  print("v_ea=list(ref_vecx)+list(ref_vecy)+list(ref_vec)")
  v_ea=np.array(v_ea)
  v_ea=v_ea.reshape(-1,3)
  print("the v_ea is "+str(v_ea))
  v_eb=[[1,0,0],[0,1,0],[0,0,1]]
  v_eb=np.array(v_eb)
  print("the v_eb is "+str(v_eb))
  print("shift v_ea to v_eb")
  cor=np.mat(cor)
  #print(cor)
  v_ea=np.mat(v_ea)
  v_eb=np.mat(v_eb)
  rot_cor=cor*np.linalg.inv(v_ea)*v_eb
  M_rot=np.linalg.inv(v_ea)*v_eb
  print("\nstep3: plus back the cor_of_c_atm_id")
  goal_cor=rot_cor+cor_of_c_atm_id
  print("reset the cor of atoms")
  #print(goal_cor)
  atoms.set_positions(goal_cor)
  return atoms,M_rot
#-------------------------------------------------------------------------

def vect_for_most_sparse_locally(atoms,c_atm_id):
  R=3
  M_V=[] 
  for i in range(0,72):
    for j in range(-18,18):
      thita=2*np.pi/72*i
      psi=np.pi/36*j
      dx=R*np.sin(psi)*np.cos(thita)
      dy=R*np.sin(psi)*np.sin(thita)
      dz=R*np.cos(psi)
      x=atoms.get_positions()[c_atm_id][0]+dx
      y=atoms.get_positions()[c_atm_id][1]+dy
      z=atoms.get_positions()[c_atm_id][2]+dz
      cor_test=[x,y,z]
      cor_test=np.array(cor_test)
      V=0
      for atm_id in range(len(atoms)):
        dis=np.linalg.norm(cor_test-atoms.get_positions()[atm_id]) 
        V=V+1/dis
      M_V=np.append(M_V,[x,y,z,V])
  M_V=M_V.reshape(-1,4)
  M_V=M_V[np.argsort(M_V[:,3])]
  print(M_V)
  cor_choose=M_V[0,0:3]

  #V_dis=[]
  #for dis in range(1,5,0.1):
  #  V=0
  #  pos=dis*vect+atoms.get_positions(c_atm_id)
  #  for atm_id in range(len(atoms)):
  #    r=np.linalg.norm(pos-toms.get_distance(atm_id))
  #    Vlj=np.power(2/atoms.get_distance(atm_id,c_atm_id,mic=True)
  #    V=V+Vlj
  vect=cor_choose-atoms.get_positions()[c_atm_id]
  vect=vect/np.linalg.norm(vect) 
  return vect


#---------------------for adb----------------------------------------


def vect_for_most_sparse_locally_v2(atoms,c_atm_id):
  
  cors=atoms.get_positions()
  cor0=atoms.get_positions()[c_atm_id]
  

  D,Dlen=get_distances(cor0,cors)
  D=D[0]
  Dlen=Dlen[0]
 
  char_vec=0
  for i in range(len(Dlen)):
    if Dlen[i]!=0:
      vec=D[i]/Dlen[i]
      char_vec = char_vec + vec
  
  vect=-char_vec

  vect=vect/np.linalg.norm(vect)
  return vect


def chaosid_2_orderid(chaosid):
  ids_content=chaosid.split(" ")
  orderid=[]
  for id_content in ids_content:
    if '-' in id_content:
      id_beg=id_content.split('-')[0]
      id_end=id_content.split('-')[1]
      ids=list(range(int(id_beg),int(id_end)+1))
    else:
      ids=int(id_content)
    orderid=np.append(orderid,ids)
  #orderid=orderid.astype(int)  
  orderid=np.unique(orderid)
  orderid=np.sort(orderid)
  return orderid

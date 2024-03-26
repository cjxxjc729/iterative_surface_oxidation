#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11


import sys
from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry
import numpy as np
from mimic_functions import *

from file_operation import *
from switch_adb import *
from structure_opt import *
from ovito_operation import *
from opt_by_orr_jobname_opt import *
from opt_by_orr_jobname_cal_dEv_and_TDI import *

import os
import time
import re
import random
import json


def indentify_adb_role(catalytic_name, atm_ids_adb, atoms):

  '''
  判断cluster中的吸附质是什么
  ref adb names are wrriten in cat_rea_adbs
  atm_ids_adb is read from *adsorbates.txt
  20240305本次修改仍然是权宜之计。后面要找到更general的分子识别方法 
  '''

  #1 ini
  ref_adb_names = resolve_cat_rea_adbs(catalytic_name)
  #print("ref_adb_terms=",ref_adb_terms)
  atoms_adb = atoms[atm_ids_adb]
  elements=np.sort(atoms_adb.get_chemical_symbols())
  adb_name=''.join(elements)

  #2 filter1
  roles=[]
  for ref_adb_name in ref_adb_names:
        
    #print("ref_adb_name=",ref_adb_name)
    ref_adb_elements=np.sort(list(ref_adb_name))
    #print("elements=",elements)

    if np.array_equal(elements, ref_adb_elements):
      #print('find the ref_adb_name, which is', ref_adb_name)
      roles.append(ref_adb_name)
      break

  #2.2 filter2
  #此处简单的设定, 所有只含有一个元素的cluster 都可以作为blank. 这很正常, 因为有时候, 比如O既可以作为吸附位, 也可以作为吸附质.
  #if len(elements)==1:
  #效果不好 改成
  if len(elements)==1 and len(roles)==0:
    roles.append('blank')

  return roles


def resolve_cat_rea_adbs(catalytic_reaction_name):

  '''
  解析文件中的吸附质
  '''

  #ref_adb_names=np.loadtxt(f_cat_rea_adbs,dtype='str')
  with open(f_reaction_template, 'r') as file:
    data = json.load(file)

  # 搜索名为catalytic_reaction_name的字典
  cat_rea_dict = data[catalytic_reaction_name]
#next((item for item in data if item["name"] == catalytic_reaction_name), None)

  # 检查是否找到了cat_rea_dict字典并打印结果
  if cat_rea_dict:
    ref_adb_terms = cat_rea_dict['rea_path']
    ref_adb_names = [] 
    
    for term in ref_adb_terms:
      ref_adb_name = [x for x in term if '*' in x][0].split("*")[0]
      if len(ref_adb_name)==0:
        ref_adb_name='blank'
      ref_adb_names.append(ref_adb_name)
      
    ref_adb_names = np.unique(ref_adb_names)
  else:
    print("未找到名为"+catalytic_reaction_name+"的字典。") 
    ref_adb_names = []

  return ref_adb_names


def fill_in_adsorbates(role, catalytic_reaction_name, atoms, atm_ids_adb):

  '''
  在已知role(adb0)的情况, 根据catalytic_reaction_name补齐剩下的吸附质
  '''
  input_dic = parse_input_file('input.parameters')
  dir_4_ref_str = input_dic['fid_ref_str']

  #1 根据role获all adbs
  adb0 = role
  adb_names = resolve_cat_rea_adbs(catalytic_reaction_name)
      
  #2
  atm_ids_adb_of_atoms_adb0 = atm_ids_adb
  adb0_name = adb0 
  atoms_adb0 = atoms  

  dict_atoms_by_adb_name = {}

  print("adb_names=",adb_names)
  for adb_name in adb_names:
    #通过switch_adb函数将adb0转成adb，同时建立起现有atoms0_adb和原有atoms_adb0的链接
    #注意当role为blank 或者不是blank时 有所不同
    if role != 'blank':    
      atoms_adb = switch_adb(adb0_name,adb_name,atoms_adb0,atm_ids_adb_of_atoms_adb0,dir_4_ref_str)
    else:
      atoms_adb = add_adb(             adb_name,atoms_adb0,atm_ids_adb_of_atoms_adb0,dir_4_ref_str) 
  
    dict_atoms_by_adb_name[adb_name] = atoms_adb
    

    print("***", adb_name,"***")
    print("atoms_adb = ", atoms_adb)

  
  return dict_atoms_by_adb_name, adb_names


def make_core(atoms, atm_ids_adb, adb_name, nl):

  #1 get c_atm_id
  indices, offset = nl.get_neighbors(atm_ids_adb[0])
  atm_ids_in_core = sort_the_indices(atoms, atm_ids_adb[0], indices)[:,0]
  atm_ids_in_core = [int(x) for x in atm_ids_in_core]
 
  #print("1: atm_ids_adb =",atm_ids_adb)
  #print("1: atm_ids_in_core=", atm_ids_in_core)

  if adb_name == 'blank':
    c_atm_id = atm_ids_adb[0]
  else:
    for atm_id in atm_ids_in_core:
      if atm_id not in atm_ids_adb:
        c_atm_id = atm_id
        break    

  #2 get atoms_core and atm_ids_in_core
  #print("atm_ids_in_core=",atm_ids_in_core)
  #print("c_atm_id=",c_atm_id)
  atm_ids_in_core = [x for x in atm_ids_in_core if x != c_atm_id] # remove(c_atm_id)
  atm_ids_in_core = [x for x in atm_ids_in_core if x not in atm_ids_adb]  
  atm_ids_in_core, true_indices = np.unique(atm_ids_in_core, return_index=True)
  atm_ids_in_core = atm_ids_in_core[np.argsort(true_indices)]

  atm_ids_left_in_core = atm_ids_in_core.copy()
  print("atm_ids_left_in_core=",atm_ids_left_in_core)
  print("c_atm_id = ",c_atm_id)  

  atoms_core = atoms[[c_atm_id]] + atoms[atm_ids_left_in_core] #+ atoms[atm_ids_adb] = atoms_with_adb
 
  atm_ids_in_core = [c_atm_id] + list(atm_ids_left_in_core) + list(atm_ids_adb)


  return atoms_core, atm_ids_in_core, c_atm_id


def site_selection(atoms, data_atmid_clusterid, catalytic_reaction_name, cutoff=8):

  #ini
  nl = ini_nl(atoms,set_cutoff=cutoff)
  uniq_cluster_id_sel = []
  uniq_clusterid = np.unique(data_atmid_clusterid[:,1])  #只选取uniq cluster id， 按照 cluster来选， 最终的结果也是cluster id

  count = -1
  while len(data_atmid_clusterid) != 0:

    count += 1

    #print("---time =", count)
    uniq_clusterid = np.unique(data_atmid_clusterid[:,1])
    #print("length of uniq_clusterid is ", len(uniq_clusterid)) 


    #1  判断：有特定的role才能选取。如果没有role， 就把对应的cluster删去, and append it into index_4_bad_cluster
    roles = []
    index_4_bad_cluster=[]

    while len(roles) == 0 and len(uniq_clusterid) != 0 :

      #1.1      
      cluster_id_sel = random.choice(uniq_clusterid)
      idx = np.where(data_atmid_clusterid[:,1]==cluster_id_sel)[0]
      atm_ids_adb_4_reader = data_atmid_clusterid[:,0][idx]
      atm_ids_adb = np.array([int(x-1) for x in atm_ids_adb_4_reader])

      #1.2
      roles = indentify_adb_role(catalytic_reaction_name, atm_ids_adb, atoms)
      print("roles=",roles)

      #1.3
      if len(roles) == 0:
        index_4_bad_cluster = np.append(index_4_bad_cluster,atm_ids_adb_4_reader)
        uniq_clusterid = list(uniq_clusterid)
        uniq_clusterid.remove(cluster_id_sel)
        uniq_clusterid = np.array(uniq_clusterid)


    uniq_cluster_id_sel.append(int(cluster_id_sel))

    represent_site_id = atm_ids_adb[0]
    indices, offsets = nl.get_neighbors(represent_site_id)

    #print("indices too close", indices+1)
    #print("indices not cluter", list(index_4_bad_cluster))
    indices_to_filter = list(index_4_bad_cluster) + list(indices+1)
    #print("indices to filter",indices_to_filter) 

    # Use a boolean mask to filter out rows where the first column matches any element in the list
    mask = ~np.isin(data_atmid_clusterid[:, 0], indices_to_filter)

    # Apply the mask to the array
    data_atmid_clusterid = data_atmid_clusterid[mask]

  return uniq_cluster_id_sel


def find_f_str_of_last_step():

  # 设定目录路径
  directory_path = './relax_trj/'

  # 正则表达式匹配形如“纯数字.xyz”的文件名
  pattern = r'^\d+\.xyz$'

  # 列出目录中所有符合条件的文件
  matching_files = [f for f in os.listdir(directory_path) if re.match(pattern, f)]
  matching_files.sort()
 
  f_str_of_last_step = matching_files[-1]
  step_index_str = f_str_of_last_step.split('.xyz')[0].split('/')[-1]
  step_index = int(step_index_str)

  return f_str_of_last_step, step_index


def delete_files_starting_with_final():

  folder_path='./relax_trj/'   

  # 检查文件夹路径是否存在
  if not os.path.isdir(folder_path):
    print(f"The folder {folder_path} does not exist.")
    return
    
  # 设置一个标志，跟踪是否找到并删除了文件
  file_deleted = False
    
  # 遍历文件夹中的所有文件
  for filename in os.listdir(folder_path):
    # 检查文件名是否以'final'开头
    if filename.startswith('final'):
    # 构造文件的完整路径
      file_path = os.path.join(folder_path, filename)
      # 删除文件
      os.remove(file_path)
      print(f"Deleted file: {file_path}")
      file_deleted = True
    
  # 如果没有找到以'final'开头的文件，则输出信息
  if not file_deleted:
    print("No files starting with 'final' were found to delete.")



if __name__ == "__main__":
   
  input_dic = parse_input_file('input.parameters')

  f_ZPE = input_dic['f_ZPE']
  f_Eaq = input_dic['f_Eaq']
  f_reaction_template = input_dic['f_reaction_template']   
  
  restart_signal = input_dic['restart'] 
  cal_tot_signal = input_dic['cal_tot']

  #---------read input ------------#
  ocutoff = float(input_dic['ocutoff'])
  icutoff = float(input_dic['icutoff'])
  f_str = input_dic['f_str']
  fid_ref_str = input_dic['fid_ref_str'] 

  #f_cat_rea_adbs='./ref_str/cat_rea_adbs.txt'  #cat_rea_adbs.txt文件记录了催化中间体
  catalytic_reaction_name = input_dic['catalytic_reaction_name']
  #--------------------------------#
  
  #--------input extended----------#
  if cal_tot_signal == 'False':
    nstep = int(input_dic['nstep'])
  
    if restart_signal == 'False':

      mkdir('relax_trj')
      f_str0 = './relax_trj/0000.xyz'  #starts from 0000.xyz
      shutil.copy(f_str, f_str0)   #将输入的f_xyz 拷贝到目录./relax_trj/中,并设计其为xyz文件
      step_beg = 0
      step_end = step_beg+nstep
    else:
      print("restat mode") 
      #找到最后一步
      delete_files_starting_with_final()  #删除'./relax_trj/'文件夹中 final 打头的文件
      f_str_of_last_step, step_index = find_f_str_of_last_step()
      restart_beg_index = str(step_index+1).zfill(4)
      print("f_str_of_last_step, step_index", f_str_of_last_step, step_index)    
      f_str0 = './relax_trj/'+restart_beg_index+'.xyz'
      shutil.copy('./relax_trj/'+f_str_of_last_step, f_str0)   #将输入的f_xyz 拷贝到目录./relax_trj/中,并设计其为xyz文件
      step_beg = step_index+1
      step_end = step_beg+nstep
  else:
    delete_files_starting_with_final()  #删除'./relax_trj/'文件夹中 final 打头的文件
    f_str_of_last_step, step_index = find_f_str_of_last_step()
    f_str0 = './relax_trj/final.xyz'
    shutil.copy('./relax_trj/'+f_str_of_last_step, f_str0) 
    step_beg = 0
    nstep = 1
    step_end = step_beg+nstep
  #--------------------------------#

  surf_ids_link={}

  for step in range(step_beg,step_end):

    formed_step = str(step).zfill(4)
    print("///////step: ",formed_step,"///////\n")
    
    '''
    1. 表面文件创建
    '''
    '''
    1.0.1创建并读取表面cluster和原子
    '''
    if cal_tot_signal == 'True':
      formed_step='final'

    f_str0 = './relax_trj/'+formed_step+'.xyz'

    clear_gas_mol(f_str0)  # clear gas mol first 

    atoms = read(f_str0)
    construct_surf_clusters(f_str0)  #通过脚本surface_cluster_reconiser.py 判断表面原子和cluster
    data_atmid_clusterid = np.loadtxt('./relax_trj/'+formed_step+'_adsorbates.txt').astype(int)  #读取文件${formed_step}_adsorbates.txt， 第一列是atm_id，第二列是clusterid. 通过astype(int)将其变成int

    '''
    1.0.2其他参数
    '''
    nl = ini_nl(atoms,set_cutoff=ocutoff)  #根据ocutoff计算nl
    #judge whether to calcualte tot TOF
    if cal_tot_signal == 'False':
      uniq_cluster_id_sel =  site_selection(atoms, data_atmid_clusterid, catalytic_reaction_name , cutoff=icutoff*3) #选取cluster, 要求距离10以上，取尽可能多的表面位
      fid_str_opt_center = 'str_opt_center'
    else:
      uniq_cluster_id_sel = np.unique(data_atmid_clusterid[:,1])
      fid_str_opt_center = 'str_opt_center.'+str(step_index).zfill(4)

    np.savetxt('./relax_trj/'+formed_step+'_clusterid_sel.txt', uniq_cluster_id_sel, fmt='%d') #save the selected clusterid
    atm_ids_in_core_tot = []   # atm_ids_in_core_tot 为所有的cutmodel里面的atm_ids, 以便后期atom_shell
    mkdir(fid_str_opt_center) #创建str_opt_center

    '''
    遍历所有的uniq_cluster_id_sel
    '''
    for cluster_id in uniq_cluster_id_sel:  #逐一遍历uniq_cluster_id_sel 并创建cut model

      print("----cluter_id---",cluster_id,'----\n')  

      '''
      1.1. 求出cluser_id所对应的atm_ids_adb,也即吸附质的atm_ids。如果说role = blank，则这个atm_ids_adb就是吸附位的atm_id
      '''
      idx = np.where(data_atmid_clusterid[:,1]==cluster_id)[0]
      atm_ids_adb = data_atmid_clusterid[:,0][idx] - 1  #
      #atm_ids_adb_4_reader = [int(x) for x in atm_ids_adb_4_reader]
      #atm_ids_adb = np.array(atm_ids_adb_4_reader)-1
      #print("atm_ids_adb=",atm_ids_adb)
      
      '''
      1.2. 得到cluster对应的possible_adb_names, (example: possible_adb_names=['blank','O'])
      修改: 实践中考虑ORR有多个吸附质容易出问题。故而暂不考虑那么多吸附质
      '''
      possible_adb_names = indentify_adb_role(catalytic_reaction_name, atm_ids_adb, atoms) #f_cat_rea_adbs是cat_rea_adbs文件位置， 里面有反应所涉及的吸附质. indentify_adb_role通过已有的吸附质， 吸附质id(atm_ids_adb), 以及atoms来判断possible_adb_names
      #将来吸附质识别更成熟后，这里indentify_adb_role要修改
      print("possible_adb_names=",possible_adb_names)

      for adb_name0 in possible_adb_names: 
 	
        '''		  
	1.3 逐一遍历adb_names, 构建吸附位和role依赖的的 O OH OOH和blank.cif 给在fid_str_opt_center文件夹中。  
        '''
        #创建atoms_core.atoms_core是去掉吸附质的版本， atm_ids_in_core是未去掉吸附质的版本. c_atm_id是经由role求得的吸附位id. 
        atoms_core, atm_ids_in_core, c_atm_id = make_core(atoms, atm_ids_adb, adb_name0, nl)  

        #将atm_ids_in_core 加入atm_ids_in_core_tot 以便后期计算atom_shell. So it is important that atm_ids_in_core是未去掉吸附质的版本
        atm_ids_in_core_tot = np.append(atm_ids_in_core_tot, atm_ids_in_core) 

        #construct the fix_list, 因为我们是将吸附质放在末尾，且总是根据blank来做fix_list,因此各个吸附质的fix_list相同
        fix_list = gen_fix_list(atoms_core, icutoff) 

        '''
	1.3.1 添加吸附质, 由于atoms_core的特点， 吸附位位于0， 所以可以方便的写[0]. 吸附质给在atoms的尾巴. 而后分别给出各个吸附质cutmodel 的f_xyz和f_fix_list

        '''
        #创建通过吸附质名字所得到atoms的字典dict_atoms_by_adb_name . adb_names1 是f_cat_rea_adbs中读出的所有吸附质的名字
        dict_atoms_by_adb_name, adb_names1 = fill_in_adsorbates('blank', catalytic_reaction_name, atoms_core, [0]) 
        #测试，考虑原始的结构可能会合理一些。 建议保留
        if adb_name0 != 'blank':
          atoms_core_with_adb_name0 = atoms_core + atoms[atm_ids_adb]
          dict_atoms_by_adb_name[adb_name0] = atoms_core_with_adb_name0
          all_adb_names = adb_names1[(adb_names1 != 'blank') & (adb_names1 != adb_name0)]
 
          natm_of_atoms_core_with_adb_name0 = len(atoms_core_with_adb_name0)           
          natm_of_adb0 = len(atm_ids_adb)
          atm_ids_adb_of_atoms_core_with_adb_name0 = range(natm_of_atoms_core_with_adb_name0 - natm_of_adb0 , natm_of_atoms_core_with_adb_name0)   #from natm_of_atoms_core_with_adb_name0 - natm_of_adb0 (atm_id of the first adsoabte's atoms) to the natm_of_atoms_core_with_adb_name0 (the end)

          for adb_name in all_adb_names:
            #natm_of_atoms_core_with_adb_name0 = len(atoms_core_with_adb_name0)
            #natm_of_adb0 = len(atm_ids_adb)  
            print("atm_ids_adb_of_atoms_core_with_adb_name0=",atm_ids_adb_of_atoms_core_with_adb_name0)
            atoms_core_with_adb_name = switch_adb(adb_name0,adb_name,atoms_core_with_adb_name0, atm_ids_adb_of_atoms_core_with_adb_name0 ,dir_4_ref_str = fid_ref_str)
            dict_atoms_by_adb_name[adb_name] = atoms_core_with_adb_name 

 
        for adb_name1 in adb_names1:
          atoms_adb = dict_atoms_by_adb_name[adb_name1]   
          prefix='C'+str(cluster_id)+'.R'+adb_name0+'.'+ adb_name1
          write('./'+fid_str_opt_center+'/'+prefix+'.xyz', atoms_adb)
          np.savetxt('./'+fid_str_opt_center+'/'+prefix+'.fix_list',fix_list,fmt='%d')      

      '''
      1.4.summary in to a dictionary
      '''
      surf_ids_link[cluster_id]={}
      surf_ids_link[cluster_id]['吸附质角色'] = possible_adb_names
      surf_ids_link[cluster_id]['吸附质原子id'] = atm_ids_adb+1

    '''    
    1.5.calculate and output atoms_shell from atm_ids_in_core_tot
    '''

    if cal_tot_signal == 'False':
      atm_ids_in_shell = [x for x in list(range(len(atoms))) if x not in atm_ids_in_core_tot]
      atoms_shell = atoms[atm_ids_in_shell]

      write('./'+fid_str_opt_center+'/shell.xyz',atoms_shell)
    
    '''
    2. structure opt
    '''
    os.chdir('./'+fid_str_opt_center+'/')
    sub_opt_by_jobs() 
  
    '''
    3. 根据计算求FED，求TDI，并合并得到atoms_new_TDI
    '''
    atoms_new_TDI = select_opt_str() 
    os.chdir('../')

    #output
    f_str_at_new_step = './relax_trj/'+str(step+1).zfill(4)+'.xyz' 
    write(f_str_at_new_step, atoms_new_TDI)
    clear_gas_mol(f_str_at_new_step)		

    print("\n----summmay-----\n")  
    print(surf_ids_link)


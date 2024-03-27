#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11

import sys
import json
from chempy import balance_stoichiometry
from chempy import Substance

from collections import OrderedDict
import random
import os
import time
import re
import glob
import subprocess

from file_operation import *

# 以氨的燃烧为例 NH3 + O2 -> NO + H2O

def blankcing_reactions(adb_rea,adb_pro,filler_mols):
  '''
  配平反应方程式
  examples of adb_rea/adb_pro: '','O','OH','COOH'
  '''
  #
  adb_rea = [item.replace('*', 'Ar') for item in adb_rea]
  adb_pro = [item.replace('*', 'Ar') for item in adb_pro]


  react_fix=set(adb_rea)
  proct_fix=set(adb_pro)
  

  #
  max_attempts = 300
  attempts = 0

  while attempts < max_attempts:

    #
    react_set = react_fix.copy()
    proct_set = proct_fix.copy()
    
    #randomly added
    for specie in filler_mols:

      index = random.randint(0, 2)

      if index==1:
        react_set.add(specie)
      elif index == 0:
        proct_set.add(specie)

    #blancing process. Try many times because we have different way to keep fillers.
    try:
      #print("react_set=",react_set)
      #print("proct_set=",proct_set)
      reactants, products = balance_stoichiometry(react_set, proct_set)
      #set x1=1
      #print("reactants=",reactants)
      #print("products=",products)

      for key, value in reactants.items():
        #print("rea, value=",value)
        if not isinstance(value, (int, float, complex)):
          if 'x1' in str(value):
            #print("find x1 in value")
            x1 = 1
            eval_value = str(value).replace('x1',str(x1))
            new_value = eval(eval_value)
            reactants[key] = new_value
        #print("rea key, value(after)=",key,reactants[key])
      for key, value in products.items():
        #print("pro key, value=",value)
        if not isinstance(value, (int, float, complex)):
          if 'x1' in str(value):
            #print("find x1 in value")
            x1 = 1
            eval_value = str(value).replace('x1',str(x1))
            new_value = eval(eval_value)
            products[key] = new_value
        #print("rea key, value(after)=",key,products[key])


      #print("reactants=",reactants)
      #print("products=",products)
      test_values=[]
      for key, value in reactants.items(): 
        test_values.append(value)
      for key, value in products.items():
        test_values.append(value)
      if not any('x' in str(item) for item in test_values) and all(value == 1 for key, value in reactants.items() if 'Ar' in key) and all(value == 1 for key, value in products.items() if 'Ar' in key) :
        break
    
    except Exception as e:
      #print(f"配平尝试 {attempts + 1} fails：{e}")  
      attempts += 1

  if attempts == max_attempts:
    print("已达到最大尝试次数，未能成功执行。")
  else:
    print(f"配平成功。尝试{attempts}次")
    print("reactants=",reactants)
    print("products=",products)

  return reactants,products


def get_dic_Gadb(jobname, dic_ZPE):

  # Pattern to match files starting with 'C1027.Rblank' and ending with '.cif'
  pattern = './' + jobname +'*.log'

  # Using glob.glob() to find all files matching the pattern
  fs_log = glob.glob(pattern)
  #print("fs_str=",fs_log)

  #print("dic_ZPE=",dic_ZPE)

  dic_Gadb={}

  for f_log in fs_log:

    adb_name = f_log.split('.')[-2].split('_')[0]
    print("f_log = ", f_log)
    with open(f_log,'r') as f:
      E = f.readlines()[-1].split()[-2]
    E = float(E)

    G = E + dic_ZPE[adb_name][0] - dic_ZPE[adb_name][1] + dic_ZPE[adb_name][2]

    dic_Gadb[adb_name] = [G]


  return dic_Gadb


def dG_calculation(reactants,products,URHE,dic_Gadb, dic_Eaq):


  #1
  Gtot_rea = 0

  #print("reactants=" ,reactants)
  for key, value in reactants.items():
  
    if key.endswith('Ar'):
      
      adb_name = key.split('Ar')[0]
      
      if len(adb_name)==0:
        adb_name='blank'
      
      G = dic_Gadb[adb_name][0]

    elif key == 'e-':

      G = -URHE

    else:
     
      G = dic_Eaq[key][0]
    
    #print("value=",value)
    nG = G * float(value) 
    #print("G  value",G,value)
    #print("nG=",G,nG)              

    Gtot_rea = Gtot_rea + nG

  #2
  Gtot_pro = 0

  for key, value in products.items():

    if key.endswith('Ar'):

      adb_name = key.split('Ar')[0]
      if len(adb_name)==0:
        adb_name='blank'

      G = dic_Gadb[adb_name][0]

    elif key == 'e-':

      G = -URHE

    else:

      G = dic_Eaq[key][0]
          
    nG = G * float(value)

    Gtot_pro = Gtot_pro + nG

  dG = Gtot_pro - Gtot_rea


  return dG


def gen_FED(jobname, catalytic_reaction_name, dic_Gadb, dic_Eaq, URHE='default'):

  
  #init--------------------------------------------
  input_dic = parse_input_file('../input_parameters.json') 

  f_ZPE = input_dic['f_ZPE']
  f_Eaq = input_dic['f_Eaq']
  f_reaction_template = input_dic['f_reaction_template']

  #ref_adb_names=np.loadtxt(f_cat_rea_adbs,dtype='str')
  with open(f_reaction_template, 'r') as file:
    data = json.load(file)

  # 搜索名为catalytic_reaction_name的字典
  #cat_rea_dict = next((item for item in data if item["name"] == catalytic_reaction_name), None)
  cat_rea_dict = data[catalytic_reaction_name]

  if URHE == 'default':
    URHE = float(cat_rea_dict["URHE.default"])
  else:
    URHE = float(URHE)

  filler_mols = cat_rea_dict["fillin_mols"]
  #------------------------------------------------
  #
  rea_path = cat_rea_dict['rea_path']

  nlevel = len(rea_path)
  nstep = nlevel -1 


  FED = [0]
  for step_idx in range(nstep):
    
    adb_rea = rea_path[step_idx]
    adb_pro = rea_path[step_idx+1]

    reactants,products = blankcing_reactions(adb_rea,adb_pro,filler_mols)    
     
    dG = dG_calculation(reactants,products,URHE, dic_Gadb, dic_Eaq) 

    G = FED[-1] + dG
    FED.append(G)

  adb_list = []
  for key in dic_Gadb:
    adb_list.append(key)
 
  return FED, adb_list
    
    


if __name__ == "__main__":
 
  #1
  #adb_rea = 'O'
  #adb_pro = 'OOH'
  input_dic = parse_input_file('../input_parameters.json')
 
  jobname = sys.argv[1] 
  f_ZPE = input_dic['f_ZPE']
  f_Eaq = input_dic['f_Eaq']
  catalytic_reaction_name = input_dic['catalytic_reaction_name']
  URHE = input_dic['URHE'] #could be default

  #2
  with open(f_ZPE, 'r') as file:
    dic_ZPE = json.load(file)
  with open(f_Eaq, 'r') as file:
    dic_Eaq = json.load(file)

  dic_Gadb = get_dic_Gadb(jobname, dic_ZPE)

  print("dic_Gadb=",dic_Gadb)
  #print("dic_Eaq=",dic_Eaq)

  #filler_mols=['H+','e-','H2O']

  #reactants,products = blankcing_reactions(adb_rea,adb_pro,filler_mols)
  #print("reactants", reactants)
  #print("products", products)
  #dG = dG_calculation(reactants,products,URHE,dic_Gadb, dic_Eaq)

  #print('dG',dG)
        
  FED, adb_list = gen_FED(jobname, catalytic_reaction_name, dic_Gadb, dic_Eaq, URHE)  
  print(adb_list)
  print(FED)


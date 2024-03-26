#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11


import sys
import os
import time
import re
import numpy as np

home_dir = os.getcwd()
fids_to_cal_tof = sys.argv[1:]
basename = os.path.basename(home_dir)


results_step_tof=[]
for fid_to_cal_tof in fids_to_cal_tof:

  print("-------------", fid_to_cal_tof, "--------------------")

  step = fid_to_cal_tof

  os.chdir(fid_to_cal_tof)

  fs=os.listdir('./')
  fs.sort()
  fs_dEv=[]
  for f in fs:
      if f.endswith('.dEv'):
          fs_dEv.append(f)

  RT=0.025875

  TOF_tot = 0

  results_prefix_dEv_tof=[]
  for f_dEv in fs_dEv:

    prefix = f_dEv.split('.')[0]

    dEv = np.loadtxt(f_dEv)

    tof = np.exp(-dEv/RT)

    results_prefix_dEv_tof = np.append(results_prefix_dEv_tof,[prefix,dEv,tof])

    TOF_tot = TOF_tot+ tof

  results_prefix_dEv_tof = results_prefix_dEv_tof.reshape(-1,3)
  sort_idx = np.argsort(results_prefix_dEv_tof[:,1].astype(float))
  results_prefix_dEv_tof = results_prefix_dEv_tof[sort_idx]

  np.savetxt(home_dir+'/results_Cid_dEv_tof.txt',results_prefix_dEv_tof,fmt='%s')

  results_step_tof = np.append(results_step_tof, [step,"{:.5e}".format(TOF_tot)])


  print("TOF tot =",TOF_tot)
  os.chdir(home_dir)

results_step_tof = results_step_tof.reshape(-1,2)
with open ('../site_activity.txt','a') as f:
  np.savetxt(f,results_step_tof,fmt='%s')



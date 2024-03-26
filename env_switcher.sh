#!/bin/bash
home_dir=$(pwd)
signature=$0


env_hpc=/public1/home/sch0149/deepmd-kit-2.2.9/bin/python3.11
env_local=/home/cjx/deepmd-kit-2.2.9/bin/python3.11

current_env=$(head -1 main.py | awk -F '!' '{print $2}')


if [ $current_env == $env_hpc ]
then
  echo "switch from env_hpc to env_local"
  sed -i "s#${current_env}#${env_local}#g" *py
elif [ $current_env == $env_local ] 
then
  echo "switch from env_local to env_hpc"
  sed -i "s#${current_env}#${env_hpc}#g" *py
fi




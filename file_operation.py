#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11

import sys

import numpy as np
import os
import time
import re
import shutil


def mkdir(path):
  if os.path.exists(path):
    shutil.rmtree(path)
  os.makedirs(path, exist_ok=True)


def mkdir_if_not_exists(directory_path):

  if not os.path.exists(directory_path):
    os.makedirs(directory_path)
  else:
    print(f"Directory {directory_path} already exists.")


def parse_input_file(filename):
    variables = {}
    with open(filename, 'r') as file:
        for line in file:
            # 去除行尾的换行符并检查是否包含'='
            if '=' in line.strip():
                # 分割键和值
                key, value = line.split('=', 1)
                # 去除键和值周围的空格
                key = key.strip()
                value = value.strip()
                # 尝试将值转换为整数，如果失败则保持字符串格式
                try:
                    value = int(value)
                except ValueError:
                    # 尝试转换为浮点数，如果失败则保持字符串格式
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                # 将键值对添加到字典中
                variables[key] = value
    return variables





if __name__ == "__main__":

  f_input = 'input.parameters'

  variables = parse_input_file(f_input)

  print(variables)

  
  

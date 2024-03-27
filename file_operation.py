#!/home/cjx/deepmd-kit-2.2.9/bin/python3.11

import sys

import numpy as np
import os
import time
import re
import shutil
import json

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
   
    with open(filename, 'r') as file:
      json_data = json.load(file)

    convert_paths_to_absolute(json_data)

    with open(filename, 'w') as file:
      json.dump(json_data, file, indent=4)

    results = extract_values_with_parent_keys(json_data)
  
    return results


def convert_paths_to_absolute(obj):
    if isinstance(obj, dict):
        for key, value in list(obj.items()):
            if isinstance(value, dict) or isinstance(value, list):
                convert_paths_to_absolute(value)
            else:
                if isinstance(value, str) and '/' in value:
                    obj[key] = os.path.abspath(value)
    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            if isinstance(item, dict) or isinstance(item, list):
                convert_paths_to_absolute(item)
            else:
                if isinstance(item, str) and '/' in item:
                    obj[i] = os.path.abspath(item)


def extract_values_with_parent_keys(obj, parent_key='', result=None):
    if result is None:
        result = {}

    if isinstance(obj, dict):
        for key, value in obj.items():
            # 如果值是字典，继续递归，但不更新结果
            if isinstance(value, dict):
                extract_values_with_parent_keys(value, key, result)
            elif isinstance(value, list):
                # 对于列表，我们遍历每个元素，但传递当前的key作为父key
                for item in value:
                    extract_values_with_parent_keys(item, key, result)
            else:
                #相对路径转绝对路径
                #if isinstance(value, str) and '/' in value:
                #    value = os.path.abspath(value)
                # 对于非字典和列表的值，直接使用父key作为结果字典的key
                result[key] = value
    elif isinstance(obj, list):
        # 如果直接就是列表，对每个元素递归调用本函数
        for item in obj:
            extract_values_with_parent_keys(item, parent_key, result)

    return result





if __name__ == "__main__":

  f_input = 'input_parameters.json'

  variables = parse_input_file(f_input)

  print(variables)

  
  

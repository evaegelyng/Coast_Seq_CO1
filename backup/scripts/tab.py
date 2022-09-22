#!/usr/bin/python
import os
import sys
import pandas as pd
import numpy as np
from os import listdir

print(sys.argv)
FF = str(sys.argv[1])
temporal_dir = str(sys.argv[2])
input_file = str(sys.argv[3])
name_file = f'{input_file}.names'

print(temporal_dir)
x_files = listdir(temporal_dir)

names = np.ravel(np.array(pd.read_csv(name_file, skiprows=0)))
names = pd.Series(range(len(names)), index=names)

x_names = np.ravel( np.array(pd.read_csv(FF)) )
intersection = np.intersect1d(x_names,names.index)
x_index = names[intersection] #get the index from index table
L = len(x_index)
if(L>0):
    print(f'{L} sequences in file {FF}')
    print(intersection)
    i=0
    for XX in x_index:
        print(f'sequence {i} out of {L}')
        bashCommand = f"sed \"{XX}q;d\" {input_file} >> {FF}_grep"
        run_comm = os.system(bashCommand)
        i = i+1
else:
    print(f'empty file {FF}')
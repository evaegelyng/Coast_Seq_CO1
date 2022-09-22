#!/usr/bin/python
import os
import sys
import pandas as pd
import numpy as np
from os import listdir

FF = str(sys.argv[1])
input_file = str(sys.argv[2])
name_file = f'{input_file}.names'

#reading index file
names = np.ravel(np.array(pd.read_csv(name_file, skiprows=0)))
names = pd.Series(range(len(names)), index=names)

intersection =(FF in names.index)

if(intersection):
    x_index = names[FF] #get the index from index table
    bashCommand = f"sed \"{x_index}q;d\" {input_file}"
    run_comm = os.system(bashCommand)
    
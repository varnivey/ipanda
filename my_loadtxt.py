
# -*- coding: utf-8 -*-

### This file is a part of iPANDA package
### (in silico Pathway Activation Network Decomposition Analysis)
### The package is distributed under iPANDA license
###
### Copyright © 2016 Insilico Medicine Inc.
###
### USA, Johns Hopkins University, ETC B301,
### 1101 East 33rd St. Baltimore, MD 21218

import os
import numpy as np



def loadtxt(infile,usecols=None,dtype=float,unpack=True,\
        skiprows=0,delimiter=' '):
    '''Alternative loadtxt which ignore commas inside quotes'''

    with open(infile,'rU') as i_file:

        return np.loadtxt(generate_lines(i_file),usecols=usecols,dtype=dtype,\
                unpack=unpack,skiprows=skiprows,delimiter=delimiter)



def generate_lines(i_file):
    '''Line generator from file'''

    for line in i_file:
        line_peaces = line.split('\"')

        for i in range(len(line_peaces)):
            if i%2 == 1:
                line_peaces[i] = line_peaces[i].replace(',','_')

        new_line = '\"'.join(line_peaces)
        yield new_line










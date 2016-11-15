
### This file is a part of iPANDA package
### (in silico Pathway Activation Network Decomposition Analysis)
### The package is distributed under iPANDA license
### 
### Copyright © 2016 Insilico Medicine Inc.
###
### USA, Johns Hopkins University, ETC B301, 
### 1101 East 33rd St. Baltimore, MD 21218

import numpy as np
import os
from my_loadtxt import loadtxt
import warnings

class pathway():
    '''class representing pathway as adjacency matrix'''

    def __init__(self, pathway_contents,name):
        '''Read pathway file to matrix'''

        self.ARR_algorithms_dict = {\
            'original':self.get_original_ARRs\
            }

        self.gene_vect = np.array([],dtype=np.str)
        self.expand_gene_list(pathway_contents)
        self.name = name


    def expand_gene_list(self,new_gene_vect):
        '''Expand gene list using pathway content file'''

        gene_vect = list(self.gene_vect)

        new_genes = list(new_gene_vect)

        gene_vect.extend(new_genes)

        genes = np.unique(gene_vect)

        genes = np.array([gene for gene in genes if gene != ''])

        self.gene_vect = genes


    def get_ARRs(self, algorithm = 'dist_sum', args = {}):
        '''Calculate ARR koefficients using given algorithm
        List of availible algorithms:

        dist_sum - simple sum of pairwize distances between nodes

        '''

        alg_call = self.ARR_algorithms_dict[algorithm]
        return alg_call(**args)


    def clean_gene_vect(self,gene_dict):
        '''Remove genes which are not in ARR_list'''

        gene_vect = list(self.gene_vect)

        gene_array = list(gene_dict[self.name])

        genes_to_remove = []

        for gene in gene_vect:
            if gene_array.count(gene) == 0:
                genes_to_remove.append(gene)

        for gene in genes_to_remove:
            gene_vect.remove(gene)

        self.gene_vect = np.array(gene_vect,dtype=np.str)


    def get_original_ARRs(self,gene_dict,arr_dict):
        '''Get ARRs from reference files'''

        return self.extract_ARRs_from_text_files(gene_dict,arr_dict,False)


    def extract_ARRs_from_text_files(self, gene_dict, arr_dict, \
                sign_only=True):
        '''Get array of signs from list of ARRs for all pathways'''

        gene_array = gene_dict[self.name]
        fig_array  =  arr_dict[self.name]

        gene_array = [x for x in gene_array if x!='']
        fig_array  = [float(x) for x in fig_array if x!='']

        gene_vect = self.gene_vect
        fin_list = []

        if sign_only:
            sign_dict = {}
            sign_list = []

            for fig in fig_array:
                if fig == 0: sign_list.append( 0.)
                if fig >  0: sign_list.append( 1.)
                if fig <  0: sign_list.append(-1.)

            for i in range(len(gene_array)):
                sign_dict[gene_array[i]] = sign_list[i]

            for gene in gene_vect:
                if sign_dict.has_key(gene):
                    fin_list.append(sign_dict[gene])
                else:
                    fin_list.append(0.)

        else:
            val_dict = {}
            for i in range(len(gene_array)):
                val_dict[gene_array[i]] = fig_array[i]

            for gene in gene_vect:
                if val_dict.has_key(gene):
                    fin_list.append(val_dict[gene])
                else:
                    fin_list.append(0.)


        return fin_list


    def get_sign_array_from_text(self, ref_file):
        '''Get sign array from text file with ARRs'''

        sign_dict = {}

        with open(ref_file, 'r') as sign_file:
            for line in sign_file:
                if not line.startswith('name'):
                    l_key,l_value = line.split(',')
                    l_value = float(l_value)
                    if l_value == 0:
                        l_value = 0
                    else:
                        l_value = l_value/np.abs(l_value)

                    sign_dict[l_key] = l_value

        gene_vect = self.gene_vect
        return [sign_dict[gene] for gene in gene_vect]














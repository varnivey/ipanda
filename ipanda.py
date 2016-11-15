#!/usr/bin/python

### This file is a part of iPANDA package
### (in silico Pathway Activation Network Decomposition Analysis)
### The package is distributed under iPANDA license
###
### Copyright © 2016 Insilico Medicine Inc.
###
### USA, Johns Hopkins University, ETC B301,
### 1101 East 33rd St. Baltimore, MD 21218

import PAS_calc as pc
import os
import optparse


parser = optparse.OptionParser()
parser.add_option('-i', dest='expr_file', help='Input file')

parser.add_option('-o', dest='outfile', help='Outfile')

parser.add_option('-p', dest='pw_dir', help='Pathway_set')

(options, args) = parser.parse_args()

if options.pw_dir    == None:
    parser.print_help()
    exit()

if options.outfile   == None:
    parser.print_help()
    exit()

if options.expr_file == None:
    parser.print_help()
    exit()

pw_dir = options.pw_dir

gene_list_file = os.path.join(pw_dir,'gene_list.txt')

modules_file = os.path.join(pw_dir,'modules.txt')
pathway_contents = os.path.join(pw_dir,'pathway_contents.txt')
pathway_file = os.path.join(pw_dir,'pathway_list.txt')

calc_alg = 'oncofinder_like'
gene_file = os.path.join(pw_dir,'pathway_contents.txt')
arr_file = os.path.join(pw_dir,'pathway_akk_new_walker.txt')

expr_file=options.expr_file
outfile=options.outfile

a = pc.expr_data()


a.calculate_CNRs_from_expr_file(expr_file,delimiter='\t',remove_quotes=True)
a.take_genes_from_list_only(gene_list_file)
a.filter_CNRs_by_ttest_continious()
a.calc_CNRs_for_modules(modules_file)

b = pc.PAS_calculation(pathway_file,pathway_contents)

gene_dict = b.get_pw_vals_dict(gene_file)
arr_dict  = b.get_pw_vals_dict( arr_file)

arr_alg = 'original'
arr_alg_args = {'gene_dict':gene_dict,'arr_dict':arr_dict}

pas_alg_args = {'ARR_alg':arr_alg,'ARR_alg_args':arr_alg_args,
				'normalize':True}

b.write_PAS_to_file(a,outfile,calc_alg,pas_alg_args)





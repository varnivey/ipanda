
### This file is a part of iPANDA package
### (in silico Pathway Activation Network Decomposition Analysis)
### The package is distributed under iPANDA license
### 
### Copyright © 2016 Insilico Medicine Inc.
###
### USA, Johns Hopkins University, ETC B301, 
### 1101 East 33rd St. Baltimore, MD 21218

import warnings
warnings.simplefilter('ignore')

import os,sys
import pathway_parser as pp
import numpy as np
from scipy.stats import ttest_1samp as ttest_func2
from scipy.stats import ttest_ind as ttest_func


class PAS_calculation:
    '''Class represents methods for PAS calculation on a set of pathways'''

    def __init__(self,pathway_file,pathway_contents):
        '''Initialize data and pathways for calculation'''

        self.alg_dict = {\
        'oncofinder_like':self.oncofinder_like_PAS,\
        }

        try:
            names_list = list(np.loadtxt(pathway_file,usecols = (0,),dtype=np.str))
        except:
            names_list = [str(np.loadtxt(pathway_file,usecols = (0,),dtype=np.str))]

        pathway_list = []

        contents_dict = self.get_pw_vals_dict(pathway_contents)

        for name in names_list:

            pathway_list.append(pp.pathway(contents_dict[name],name))

        self.pathway_list = pathway_list


    def get_pw_vals_dict(self, pathway_vals):
        '''Read pathway contents file to a dictionary'''

        with open(pathway_vals,'r') as r_file:
            line1 = r_file.readline()
            line1 = line1.replace('\n','')
            line1 = line1.replace('\r','')
            pathway_names = line1.split('\t')

        all_vals = np.loadtxt(pathway_vals,skiprows=1,unpack=True,\
                dtype=np.str,delimiter='\t')

        pw_val_dict = {}

        if len(pathway_names) >= 1:
            for i,pw in enumerate(pathway_names):
                pw_val_dict[pw] = list(all_vals[i])
        else:
            for pw in pathway_names:
                pw_val_dict[pw] = list(all_vals)

        return pw_val_dict



    def write_PAS_to_file(self, data, output_file, \
            algorithm='oncofinder_like', args={}):
        '''Write PASes for given algorithm'''

        alg_call = self.alg_dict[algorithm]

        PAS_vals = alg_call(data,**args)

        pw_list = [pw.name for pw in self.pathway_list]

        self.write_pw_values_file(pw_list, data.t_measures,PAS_vals,output_file)


    def write_pw_values_file(self, pw_list, measures, vals, outfile):
        '''Write annotated table with values, sample names and pathway names'''

        if vals.ndim == 1:
            vals = vals.reshape((1,-1))

        with open(outfile,'w') as o_file:
            o_file.write('\"'+'\" \"'.join(measures)+'\"\n')
            for i in range(len(pw_list)):
                o_file.write('\"' + pw_list[i] + '\" ')
                vals_text = [str(val) for val in vals[i]]
                o_file.write(' '.join(vals_text) + '\n')


    def oncofinder_like_PAS(self,data,ARR_alg,ARR_alg_args,\
                normalize=False):
        '''Calculate oncofinder like PASes'''

        pathway_list = self.pathway_list
        genes,CNRs = data.genes, data.CNRs

        PAS_vals = []

        for pathway in pathway_list:
#            print pathway.name
#            pathway.clean_gene_vect(gene_list_file)
            ARRs = pathway.get_ARRs(ARR_alg,ARR_alg_args)
            ARRs = np.array(ARRs)
            gene_vect = pathway.gene_vect
            cur_CNRs = []
            for gene in gene_vect:
                if genes.count(gene) != 0:
                    cur_CNRs.append(CNRs[:,genes.index(gene)])
                else:
                    cur_CNRs.append(np.zeros((CNRs.shape[0]),dtype=np.float))
            cur_CNRs = np.array(cur_CNRs)
#            print len(ARRs)
#            print cur_CNRs.shape
            if normalize:
                ARRs_neg = ARRs[ARRs < 0]
                ARRs_pos = ARRs[ARRs > 0]
                ARRs_sum = max(np.sum(np.abs(ARRs_neg)),np.sum(ARRs_pos))
                if ARRs_sum != 0:
                    ARRs = 1000*ARRs/ARRs_sum

            PAS_vals.append(np.dot(cur_CNRs.T,ARRs))

        return np.array(PAS_vals)


class expr_data():
    '''Class for managing expression datasets'''

    def __init__(self):
        '''Constructor'''
        pass


    def calculate_CNRs_from_expr_file(self,expr_file,delimiter = '\t',
            type_line = False, perform_t_test = False, remove_quotes=False):
        '''Calculate CNRs using expression raw data'''


        with open(expr_file,'rU') as e_file:
            line1 = e_file.readline()
            line1 = line1.replace('\n','')
            line1 = line1.replace('\r','')
            measures = line1.split(delimiter)[1:]
            if remove_quotes:
                measures = [measure.replace('\"','') for measure in measures]
            ncols = len(measures)
            skiprows = 1

            line2 = e_file.readline()
            line2 = line2.replace('\"','')
            type_line = line2.startswith('Type')

            if type_line:
                line2   = line2.replace('\n','')
                line2   = line2.replace('\r','')
                m_types = line2.split(delimiter)[1:]
                skiprows = 2
                t_line = e_file.readline()
                t_line = t_line.replace('\"','')
            else:
                t_line = line2

            if t_line.split(delimiter)[0] == '1':
                sc = 1
            else:
                sc = 0

        expr_full = np.loadtxt(expr_file,skiprows=skiprows, unpack=True,\
                usecols = tuple(np.arange(sc+1,sc+ncols+1)),delimiter=delimiter,
                dtype=np.str)

        expr_full[expr_full == 'NA'] = '0.'
        expr_full = expr_full.astype(np.float)

        genes = list(np.loadtxt(expr_file,skiprows=skiprows,usecols=(sc,),\
                dtype=np.str, delimiter=delimiter))

        if remove_quotes:
            genes = [gene.replace('\"','') for gene in genes]

        if type_line:
            n_measures_inds   = [measures.index(m) for m in measures \
                    if m_types[measures.index(m)].startswith('normal')]
            t_measures_inds   = [measures.index(m) for m in measures \
                    if not m_types[measures.index(m)].startswith('normal')]
            t_measures  = [m for m in measures \
                    if not m_types[measures.index(m)].startswith('normal')]
            n_measures  = [m for m in measures \
                    if     m_types[measures.index(m)].startswith('normal')]

        else:
            n_measures_inds   = [measures.index(m) for m in measures \
                    if     m.startswith('Normal')]
            t_measures_inds   = [measures.index(m) for m in measures \
                    if not m.startswith('Normal')]
            t_measures  = [m for m in measures if not m.startswith('Normal')]
            n_measures  = [m for m in measures if     m.startswith('Normal')]


#        expr_full_log = np.log(expr_full)/np.log(10)

        genes = np.array(genes)
        bad_expr = expr_full <= 0
        good_inds = np.logical_not(bad_expr.sum(0))

        expr_full = expr_full[:,good_inds]
        genes = list(genes[good_inds])

        expr_full_log = np.log(expr_full)

        expr_n = expr_full_log[n_measures_inds]
        n_mean = np.mean(expr_n,0)

        expr_t = expr_full_log[t_measures_inds]
        CNRs = expr_t - n_mean

        self.expr_t = expr_t
        self.expr_n = expr_n
        self.genes = genes
        self.CNRs = CNRs
        self.t_measures = t_measures
        self.n_measures = n_measures


    def calc_CNRs_for_modules(self, modules_file,mod_cnr='mean'):
        '''Calculate CNRs for modules'''

        genes = self.genes
        CNRs = self.CNRs

        modules = []

        with open(modules_file,'r') as m_file:
            for line in m_file:
                line = line.replace('\n','')
                modules.append(line.split())

        new_CNRs = np.empty((CNRs.shape[0],len(modules)))

        for i,module in enumerate(modules):
            lines = []
            for gene in module:
                if genes.count(gene) > 0:
                    lines.append(CNRs[:,genes.index(gene)])
            if len(lines) != 0:
                if mod_cnr == 'mean':
                    new_CNRs[:,i] = np.mean(lines,0)
                elif mod_cnr == 'sum':
                    new_CNRs[:,i] = np.sum(lines,0)
            else:
                new_CNRs[:,i] = 0.

        module_names = [('mod' + str(x).zfill(4)) for x in range(len(modules))]
        genes = genes + module_names

        self.genes = genes

        self.CNRs = np.hstack((CNRs,new_CNRs))


    def take_genes_from_list_only(self,gene_list_file):
        '''Remove all genes except those from file given'''

        genes  = self.genes
        expr_t = self.expr_t
        expr_n = self.expr_n
        CNRs   = self.CNRs

        if type(gene_list_file) == str:
            with open(gene_list_file, 'r') as gl_file:
                gene_list = [line.replace('\n','') for line in gl_file]
        else:
            gene_list = list(gene_list_file)

        gene_inds = []

        for gene in gene_list:
            if genes.count(gene) > 0:
                gene_inds.append(genes.index(gene))

        self.CNRs    = CNRs[:,gene_inds]
        self.genes   = list(np.array(genes)[gene_inds])
        self.expr_n  = expr_n[:,gene_inds]
        self.expr_t  = expr_t[:,gene_inds]


    def filter_CNRs_by_single_ttest_continious(self, min_cf=0.0000001,max_cf=0.1):
        '''Filter CNRs by single sample t-test with smooth threshold'''

        if (len(self.n_measures) == 1) and (len(self.t_measures) == 1):
            print "Warning! Only 1 sample in each group"
            print "Statistical weights are off"
            return

        if not hasattr(self, 's_pval_list'):
            self.calculate_single_sample_ttest_pvals()

        CNRs = self.CNRs
        pvals = self.s_pval_list

        koefs = (np.cos((np.log(pvals) - np.log(min_cf))/(np.log(max_cf) - np.log(min_cf))*np.pi) + 1.)/2.

        koefs[pvals > max_cf] = 0.
        koefs[pvals < min_cf] = 1.

        self.s_koefs = koefs

        self.CNRs = CNRs*koefs



    def calculate_single_sample_ttest_pvals(self):
        '''Calculate p-values from single sample t-test'''

        expr_t = self.expr_t
        expr_n = self.expr_n

        mes_size, gene_size = expr_t.shape

        all_pvals = []

        for i in range(gene_size):
            n_sample  = expr_n[:,i]
            cur_pvals = []
            n_sample_r = randomize_samples(n_sample)

            cur_pvals = [ttest_func2(n_sample_r,expr_val)[1] for expr_val in expr_t[:,i]]

            all_pvals.append(cur_pvals)

        all_pvals = np.array(all_pvals).T
        all_pvals[all_pvals == 0] = sys.float_info.min

        self.s_pval_list = all_pvals



    def filter_CNRs_by_ttest_continious(self):
        '''Filter CNRs by group t-test without cutoff'''

        min_cf = 0.0000001
        max_cf = 0.1

        if self.CNRs.shape[0] == 1:
            self.filter_CNRs_by_single_ttest_continious(min_cf,max_cf)
            return

        if not hasattr(self, 'pval_list'):
            self.calculate_group_ttest_pvals()

        pval_list = self.pval_list

        self.filter_CNRs_common_continious(pval_list,min_cf,max_cf)



    def filter_CNRs_common_continious(self, val_list, min_cf, max_cf):
        '''Filter CNRs using list of values'''

        CNRs = self.CNRs
        koefs = []

        for i in range(len(self.genes)):

            cur_val = val_list[i]

            if cur_val < min_cf:
                koef = 1.
            elif cur_val > max_cf:
                koef = 0.
            else:
                koef = (np.cos((np.log(cur_val) - np.log(min_cf))/(np.log(max_cf) - np.log(min_cf))*np.pi) + 1.)/2.

            koefs.append(koef)
            CNRs[:,i] = CNRs[:,i]*koef

        self.koefs = koefs
        self.CNRs = CNRs


    def calculate_group_ttest_pvals(self):
        '''Calculate p-values from group t-test'''

        expr_t = self.expr_t
        expr_n = self.expr_n

        pval_list = []

        for i in range(len(self.genes)):


            t_sample = expr_t[:,i]
            n_sample = expr_n[:,i]

            t,cur_pval = ttest_func(randomize_samples(t_sample),randomize_samples(n_sample))

            pval_list.append(cur_pval)

        pval_list = np.array(pval_list)
        pval_list[pval_list == 0] = sys.float_info.min

        self.pval_list = list(pval_list)


def randomize_samples(sample):
    '''Randomize samples with 5% varience in case of equal samples'''

    n_samp = len(sample)

    if n_samp == 1 or len(np.unique(sample)) > 1:
        return sample

    mean = np.mean(sample)
    std = mean/20.
    if std == 0:
        std = 0.001

    return np.random.normal(mean,std,n_samp)



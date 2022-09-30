##########################################################INFO###############################################################################################
#run_codeml.py			
#
#  Script to calculate LRT values from log likelihood values obtained with codeml
#
#  Author: Diana Moreno 
#  v.1.0 september 2020
#
#  Description:
#

#
#		
#
#usage: python lrt_construction.py -in <data frame with log values alternative and null> -out <output_name>

#

#
################################################################################################################################################################



#----------------------------------------------------------------- SETTING ENVIRONMENT -----------------------------------------------------------------

#1)Import modules:

import argparse
import os
import os.path
import sys 
import pandas 
import numpy 


#2)Define input arguments:

LNL=pandas.read_csv('all_lnl_values', sep='\t', names=["Gene1", "Alt", "Gene2", "Null"])

LNL['Diff']= LNL['Alt'] - LNL['Null']

LNL['LRT']= LNL['Diff'] *2

del LNL['Gene2']


LNL.to_csv('LRT_results', sep='\t', index=False)


##########################################################INFO###############################################################################################
#run_codeml.py			
#
#  Setup script for generate CODEML null and alternative model ctl files.
#
#  Author: Diana Moreno 
#  v.1.0 January 2020
#
#  Description:
#

#
#		
#
#usage: python codeml_ctl.py -aln <codon_alignment.phy> -tree <gene_tree.newick>

#

#
################################################################################################################################################################



#----------------------------------------------------------------- SETTING ENVIRONMENT -----------------------------------------------------------------

#1)Import modules:

import argparse
import os
import os.path
import subprocess
from subprocess import check_output
import sys 
import time


#2)Define input arguments:

def get_args():
	parser = argparse.ArgumentParser(description='Input needed to build codeml ctl files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-aln', '--alignment', type=str, help='CDS alignment file in phylip format', required=True)
	parser.add_argument('-tr', '--tree', type=str, help='Gene phylogeny in newick format', required=True)

	args = parser.parse_args()
	ALIGN = args.alignment 
	TREE = args.tree
	return ALIGN, TREE

ALIGN, TREE = get_args()
	
PX=ALIGN.split("_")[0]
WORKDIR= os.getcwd()
#3)Create directories for each task

os.mkdir('null_model') 
os.mkdir('alternative_model') 


##Create codeml ctl file for NULL model 

OUTNULL=PX + '_null_model'

NULLCTL=('codeml_null.ctl')

with open (NULLCTL, 'w') as f:
	f.write ('seqfile  = ' + WORKDIR + '/' + ALIGN + '    * sequence data file name' + '\n')
	f.write ('treefile = ' +  WORKDIR + '/' + TREE + '     * tree structure file name' + '\n') 
	f.write ('outfile  = ' + OUTNULL + '  * main result file name' + '\n') 
	f.write ('\n')
	f.write ('    noisy = 9     * 0,1,2,3,9: how much rubbish on the screen' + '\n') 
	f.write ('  verbose = 1     * 1: detailed output, 0: concise output' + '\n')
	f.write ('  runmode = 0     * 0: user tree;  1: semi-automatic;  2: automatic' + '\n') 
	f.write ('                  * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise' + '\n') 
	f.write ('  seqtype = 1     * 1:codons; 2:AAs; 3:codons-->AAs' + '\n') 
	f.write ('CodonFreq = 2     * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table' + '\n') 
	f.write ('    clock = 0     * 0: no clock, unrooted tree, 1: clock, rooted tree' + '\n') 
	f.write ('   aaDist = 0     * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}' + '\n') 
	f.write ('    model = 2     * models for codons:' + '\n') #2 is used for branch-site model
	f.write ('                  * 0:one, 1:b, 2:2 or more dN/dS ratios for branches' + '\n')
	f.write ('  NSsites = 2     * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;' + '\n')#Set 2 to allows 3 categories for sites:(purifying, neutral and positive)
	f.write ('                  * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal' + '\n')
	f.write ('    icode = 0     * 0:standard genetic code; 1:mammalian mt; 2-10:see below' + '\n')
	f.write ('    Mgene = 0     * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all' + '\n')
	f.write ('\n')
	f.write ('fix_kappa = 0      * 1: kappa fixed, 0: kappa to be estimated' + '\n')
	f.write ('    kappa = 2      * initial or fixed kappa' + '\n')
	f.write ('fix_omega = 1      * 1: omega or omega_1 fixed, 0: estimate' + '\n')
	f.write ('    omega = 1      * initial or fixed omega, for codons or codon-based AAs' + '\n')
	f.write ('\n')
	f.write ('       getSE = 0          * 0: dont want them, 1: want S.E.s of estimates' + '\n') 
	f.write ('RateAncestor = 0          * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)' + '\n')  
	f.write ('  Small_Diff = .45e-6     * Default value.' + '\n')  
	f.write ('   cleandata = 1          * remove sites with ambiguity data (1:yes, 0:no)?' + '\n') 
	f.write (' fix_blength = 0          * 0: ignore, -1: random, 1: initial, 2: fixed')
	

##Create codeml ctl file for alternative model (under selection)
OUTALT=PX + '_alternative_model'
ALTCTL=('codeml_alt.ctl')

with open (ALTCTL, 'w') as f:
	f.write ('seqfile  = ' +  WORKDIR + '/' + ALIGN + '  * sequence data file name' + '\n')
	f.write ('treefile = ' +  WORKDIR + '/' + TREE + '  * tree structure file name' + '\n') 
	f.write ('outfile  = ' + OUTALT + '  * main result file name' + '\n') 
	f.write ('\n')
	f.write ('    noisy = 9       * 0,1,2,3,9: how much rubbish on the screen' + '\n') 
	f.write ('  verbose = 1       * 1: detailed output, 0: concise output' + '\n')
	f.write ('  runmode = 0       * 0: user tree;  1: semi-automatic;  2: automatic' + '\n') 
	f.write ('                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise' + '\n') 
	f.write ('  seqtype = 1       * 1:codons; 2:AAs; 3:codons-->AAs' + '\n') 
	f.write ('CodonFreq = 2       * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table' + '\n') 
	f.write ('    clock = 0       * 0: no clock, unrooted tree, 1: clock, rooted tree' + '\n') 
	f.write ('   aaDist = 0       * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}' + '\n') 
	f.write ('    model = 2       * models for codons:' + '\n') #2 is used for branch-site model
	f.write ('                    * 0:one, 1:b, 2:2 or more dN/dS ratios for branches' + '\n')
	f.write ('  NSsites = 2       * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;' + '\n')#Set 2 to allows 3 categories for sites:(purifying, neutral and positive)
	f.write ('                    * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal' + '\n')
	f.write ('    icode = 0       * 0:standard genetic code; 1:mammalian mt; 2-10:see below' + '\n')
	f.write ('    Mgene = 0       * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all' + '\n')
	f.write ('\n')
	f.write ('fix_kappa = 0       * 1: kappa fixed, 0: kappa to be estimated' + '\n')
	f.write ('    kappa = 2       * initial or fixed kappa' + '\n')
	f.write ('fix_omega = 0       * 1: omega or omega_1 fixed, 0: estimate' + '\n')
	f.write ('    omega = 1       * initial or fixed omega, for codons or codon-based AAs' + '\n')
	f.write ('\n')
	f.write ('       getSE = 0         * 0: dont want them, 1: want S.E.s of estimates' + '\n') 
	f.write ('RateAncestor = 0         * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)' + '\n')  
	f.write ('   Small_Diff = .45e-6   * Default value.' + '\n')  
	f.write ('    cleandata = 1        * remove sites with ambiguity data (1:yes, 0:no)?' + '\n') 
	f.write ('  fix_blength = 0        * 0: ignore, -1: random, 1: initial, 2: fixed')
	


subprocess.check_call(['mv', NULLCTL, 'null_model/']) 

subprocess.check_call(['mv', ALTCTL, 'alternative_model/']) 


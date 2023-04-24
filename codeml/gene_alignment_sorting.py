import os
import sys
import itertools
from collections import OrderedDict
from string import ascii_lowercase
from pathlib import Path
import shutil

#### INPUTS ####
WORK_DIR = '/lustre/scratch/npaulat/selection/alignment_groups/'
Path(WORK_DIR).mkdir(exist_ok=True)

#ALIGNMENT_LIST_FILE = sys.argv[1]
#ALIGNMENT_LIST_FILE = '/lustre/scratch/npaulat/selection/grouping_test_list'
ALIGNMENT_LIST_FILE = '/lustre/scratch/npaulat/selection/gene_alignment_file_list'

with open(ALIGNMENT_LIST_FILE) as f:
	ALIGNMENT_FILE_LIST = list(line for line in f.read().splitlines() if line)
	ALIGNMENT_FILE_LIST.sort()


#### SPECIES LISTS ####
LAURASIATHERIA = ['hg38', 'bosTau9', 'canFam4', 'equCab3', 'felCat9', 'mm10', 'susScr11']

PAIRS = ['HLmyoMyo6', 'HLmyoNig1', 'HLdipEca1', 'HLdesRot8A', 'HLmolMol2', 'HLtadBra2', 'HLsacBil1', 'HLartJam3']

RHINOLOPHIDAE = ['HLrhiAff2', 'HLrhiFer5', 'HLrhiLuc4', 'HLrhiSed2', 'HLrhiTri2']
HIPPOSIDERIDAE = ['HLaseSto2', 'HLhipCyc2', 'HLhipLar2']
PTEROPODIDAE = ['HLrouAeg4']

## excluded family lists ##
#RHINOPOMIDAE = ['HLrhiMic1'] # Collab group decided not needed, no LQ
#SACCOPTERIDAE = ['HLsacBil1', 'HLsacLep1'] # unnecessary since 1 sp in PAIRS list
#PHYLLOSTOMIDAE = ['HLartJam3', 'HLdesRot8A', 'HLdipEca1', 'HLphyDis3', 'HLphyHas1'] # unnecessary since 3 sp are in the PAIRS list
#MOLOSSIDAE = ['HLmolMol2', 'HLtadBra2'] # this is redundant since both are in the PAIRS list
#VESPERTILIONIDAE = ['HLantPal2', 'HLmyoMyo6', 'HLmyoNig1', 'HLpipKuh2'] # unnecessary since 2 sp in PAIRS list

## list of clade lists ##
#CLADE_CHECKLISTS = [RHINOLOPHIDAE, HIPPOSIDERIDAE, RHINOPOMIDAE, PTEROPODIDAE, SACCOPTERIDAE, PHYLLOSTOMIDAE, MOLOSSIDAE, VESPERTILIONIDAE, LAURASIATHERIA]
CLADE_CHECKLISTS = [RHINOLOPHIDAE, HIPPOSIDERIDAE, PTEROPODIDAE, LAURASIATHERIA]


#### CHECKLIST FUNCTIONS ####
def all_pairs(ALIGNMENT_SP, PAIR_LIST):
	SET_A = set(ALIGNMENT_SP)
	SET_B = set(PAIR_LIST)
	if SET_B.issubset(SET_A):
		return True
	else:
		return False

def enough_clades(ALIGNMENT_SP, CLADE_LISTS):
	FILE_SET = set(ALIGNMENT_SP)
	CHECK = 'No'
	for FAMILY_LIST in CLADE_LISTS:
		CLADE_SET = set(FAMILY_LIST)
		if bool(FILE_SET & CLADE_SET):
			CHECK = 'Yes'
			continue
		else:
			CHECK = 'No'
			break
	return CHECK

def get_headers(ALIGNMENT_FILE):
	SPECIES_LIST = []
	with open(ALIGNMENT_FILE, 'r') as fd:
		for line in fd.readlines():
			if '>' in line:
				NAME = line.strip()[1:] # to remove the '>'
				SPECIES_LIST.append(NAME)
	return SPECIES_LIST


#### RUN CHECKLIST FUNCTIONS ####
# 1. Initialize objects:
#	a) List of alignment files with all 30 species
#	b) List of alignment files <30 species that have required minimum species
#	c) Dictionary for alignment files (b) with their set of species
# 2. Loop through alignment files and sort based on checklist requirements
# 3. If an alignment meets the minimum requirements:
#	a) Append to alignment list
#	b) Append alignment file name to dictionary and set species list (delimited text) as value

COMPLETE_ALIGNMENT_LIST = []
GOOD_ALIGNMENT_LIST = []
FILE_SETS =  {}
for ALIGN_FILE in ALIGNMENT_FILE_LIST:
	AL_SP = get_headers(ALIGN_FILE)
	AL_SP.sort()
	SP_COUNT = len(AL_SP)
	if SP_COUNT == 30:
		#FILENAME = os.path.basename(FILE)
		COMPLETE_ALIGNMENT_LIST.append(ALIGN_FILE)
	else:
		if SP_COUNT >= 12:
			if all_pairs(AL_SP, PAIRS):
				ENOUGH_CLADES = enough_clades(AL_SP, CLADE_CHECKLISTS)
				if ENOUGH_CLADES == 'Yes':
					GOOD_ALIGNMENT_LIST.append(ALIGN_FILE)
					FILE_SP = ','.join(AL_SP)
					FILE_SETS[ALIGN_FILE] = FILE_SP
###

#put filename and headers (species) into sets and group same sets
#FILE_SETS =  {}
#FILE_SP = ','.join(AL_SP)
#FILE_SETS[FILE] = FILE_SP

## Invert dictionary so that the species list object is the key, and the associated value is a list of gene alignment files with that specific species list ##
INV_SETS = {}
for k, v in FILE_SETS.items():
	INV_SETS[v] = INV_SETS.get(v, []) + [k]

## or try
# from collections import defaultdict
# SP_SETS = defaultdict(list)
# for k, v in FILE_SETS.items():
	# SP_SETS[v].append(k)

## Sort inverted dictionary by key value (largest species list object to smallest) AND alphabetically sort the associated values (each key's associated list of gene alignment files) ##
temp_dict = {key : sorted(INV_SETS[key]) for key in sorted(INV_SETS, reverse=True)}
SP_SETS = OrderedDict(temp_dict)

## Examples in comments ##
#d = {'lion,elk': ['ECHO', 'DELTA'], 'lion,elk,fruit': ['ECHO', 'DELTA'], 'lion': ['ECHO', 'DELTA']}
#d1 = {key : sorted(d[key]) for key in sorted(d, reverse=True)}
#	{'lion': ['ECHO', 'DELTA'], 'lion,elk': ['ECHO', 'DELTA'], 'lion,elk,fruit': ['ECHO', 'DELTA']}
#OrderedDict(d1)
#	OrderedDict([('lion,elk,fruit', ['ECHO', 'DELTA']), ('lion,elk', ['ECHO', 'DELTA']), ('lion', ['ECHO', 'DELTA'])])
##

# To just order keys by largest species list
#SP_SETS = sorted(INV_SETS.items(), reverse=True)


#### Use inverted dictionary to create subdirectories for unique species lists and copy gene alignments to the associated subdirectory ####
TREE_COUNT = len(SP_SETS)

def iter_alpha_string():
	for i in itertools.count(1):
		for s in itertools.product(ascii_lowercase, repeat=i):
			yield ''.join(s)

FOLDER_STRINGS = []
for s in itertools.islice(iter_alpha_string(), TREE_COUNT):
	FOLDER_STRINGS.append(s)

FOLDER_TREES = []
TREES = list(SP_SETS)
for SUB, TREE in zip(FOLDER_STRINGS, TREES):
	NUM_SP = len(TREE.split(','))
	TREE_INFO = SUB + '\t' + str(NUM_SP) + '\t' + TREE
	FOLDER_TREES.append(TREE_INFO)
	SUBDIR = os.path.join(WORK_DIR, SUB)
	Path(SUBDIR).mkdir(exist_ok=True)
	FILES = list(SP_SETS[TREE])
	for FILE in FILES:
		BASE_FILE = os.path.basename(FILE)
		CP_FILE = os.path.join(SUBDIR, BASE_FILE)
		shutil.copy2(FILE, CP_FILE)
		FOLDER_TREES.append(CP_FILE)

FOLDER_TREES.append('og_tree\t30\tall_species')
SUBDIR = os.path.join(WORK_DIR, 'og_tree')
Path(SUBDIR).mkdir(exist_ok=True)
for FILE in COMPLETE_ALIGNMENT_LIST:
	BASE_FILE = os.path.basename(FILE)
	CP_FILE = os.path.join(SUBDIR, BASE_FILE)
	shutil.copy2(FILE, CP_FILE)
	FOLDER_TREES.append(CP_FILE)

## Output alignment subdirectory and grouping info ##
OUTPUT_LIST_FILE = os.path.join(WORK_DIR, 'alignment_tree_grouping.txt')
with open(OUTPUT_LIST_FILE, 'w') as f:
	for line in FOLDER_TREES:
		f.write("%s\n" % line)

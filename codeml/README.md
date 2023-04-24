# Running CodeML

Workflow for creating CodeML input files and running CodeML.

Originally I was supposed to run CodeML just for genes in the DNA repair and TE defense gene lists, but there were a couple of issues with that:
1. The *a priori* lists are based on human/mouse gene function annotations, which are incomplete and biased, since we are looking at genes in bats. So these lists will not contain all relevant genes for DNA repair/TE defense/immunity/etm, and genes that have relevant functions in human/mouse might not in bats.
2. Of the 17k+ genes, not all gene alignments contain sequences from all 30 species, and CodeML requires input trees with only species that are in the sequence file, so of several hundred genes on the two lists, less than half had all 30 species. Since this was a problem that affected ~6k genes, I had to make a dynamic R script to make appropriately rooted and labeled trees for each gene alignment * each hypothesis.

## Setting up

### Sorting Gene Alignments
To identify gene alignment files with <30 species represented, and exclude any not meeting the basic species representation criteria, I made the *gene_alignment_sorting.py* script. It produces the output alignment_tree_grouping.txt, which is a list of alignment files that meet filtering standards and are grouped by the species present in the alignment. This is mainly informational, the main point is to get the filtered list of alignment files.
  * Minimum species criteria:
    * All species used in the longevity comparison pairs (8 total)
    * At least one representative from each of the 8 bat families included in the study (excluding Rhinopomidae, due to lack of relevant longevity info)
    * At least one outgroup Laurasiatherian species (of 7 total)

Note: alignments must have at least 12 species total to possibly meet the above criteria.

### Input Trees
CodeML requires a rooted Newick tree with labeled nodes/branches for testing differing node/branch selection models/rates (see [PAML/codeml documentation](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf) pg. 15 *branch or node labels* section).
  * Null model is that all taxa have same selection signal (foreground = background; all nodes use the default #0 w rate).
  * Alternative model is that different branches have different selection signals (foreground =/= background; at least one node/branch uses alternate #1 rate label or $1 node label).
  * To test different selection rate hypotheses, you need differently labeled Newick tree files.

**NOTE:** CodeML requires that **ONLY** the species present in the gene alignment are in the input tree. So if running on multiple gene alignments (yes), you will need multiple clipped trees that reflect whatever combination of species actually have sequences in the given gene alignment file.

  1. Make or obtain a phylogenetic tree with branch length information, convert to Newick format.
      * In this case, the [phylogenetic tree](https://hackmd.io/6RsnHiNGTL69qaed5xcg1g) Diana and I used was made by Graham (need to get methods info). Contains bootstrap support values, though not needed here.

  2. Reroot tree if necessary (yes).
        * To make the [**rooted_tree.nwk**](https://hackmd.io/St77zA2MSXyN21bNMUAbOw) (or to modify it) in R, use ape and ggtree.
      
            ```
            library(ape)
            library(treeio)
            library(ggtree)
            library(tidytree)

            tree=read.tree("forcedtop_tree.tre")
            tree$node.label <- NULL
            tree_info <- as_tibble(tree)
            step1 <- tree_info %>% filter(!grepl('HL|100', label))
            out_info <- step1 %>% filter(!is.na(label))
            out_nodes <- range(out_info$node)

            #bat_rooted=(root(tree, 12:18, resolve.root = TRUE))
            bat_rooted=(root(tree, c(out_nodes[1]:out_nodes[2]), resolve.root = TRUE))
            ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
            ggtree(bat_rooted) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
            is.rooted(bat_rooted) # <- programatic check, should return "TRUE"
            plot(bat_rooted, show.tip.label = TRUE)
            write.tree(bat_rooted, file="rooted_tree.nwk")
            ```
        * To make an unrooted version (**unrooted_tree.nwk**) of the above tree:
            ```
            bat2=(unroot(bat_rooted))
            #plot(bat2, show.node.label = TRUE) # <- visual check
            #is.rooted(bat2) # <- programatic check, should return "FALSE"
            write.tree(bat2, file="unrooted_tree.nwk")
            ```
  3. Make labeled trees, one for each hypothesis to test, using [**reroot_label_nwk_trees2.R**](https://hackmd.io/sJelYIlLRq6-rbZWYbFpVg).
     * To create the five different foreground node labeled trees for each gene alignment file that does not have all 30 species, run **reroot_label_nwk_trees2.R** on the command line, looping through the individual gene alignment clipped tree files (e.g. ENST00000374900.FAAH2.cleanLb_hmm_manual.fasta.nwk).
        ```
        cd /lustre/scratch/npaulat/selection/trees/
        for i in *.nwk; do Rscript /lustre/scratch/npaulat/selection/reroot_label_nwk_trees2.R ${i}; done
        ``` 
       * If all of your files contain the same species, you only need one tree for each hypothesis (5), as opposed to gene-specific trees for each hypothesis(# of genes * 5).
     * Hypotheses:
       1. Bats have different selection rate/pressure than other mammals.
           * Place #1 label on the ancestral bat node (can check in tree view to make sure it's on the correct node). 
       2. Long-lived bat species have different selection signal from other (i.e. short-lived) bat species/outgroups.
           * Move #1 from the ancestral bat node to the branch of each long-lived species (e.g. HLmyoMyo6 #1:<branch length>).
       3. Short-lived bat species have different signal from other (i.e. long-lived) bat species/outgroups.
           * Move #1 from the ancestral bat node to the branch of each short-lived species (e.g. HLmyoNig1 #1:<branch length>).
       4. Vespertilionid bats have a different selection signal than other bats/outgroups.
           * Move #1 from the ancestral bat node to the base of Vespertilionidae.
       5.  Rhinolophid bats have a different selection signal than other bats/outgroups.
           * Move #1 from the ancestral bat node to the base of Rhinolophidae. 
<br>
  4. Spot check a few Newick files to check that the trees are properly labeled.
<br>

### Other CodeML Input Files
  5. [**ctl_files.py**](https://hackmd.io/sRHc3iUpQ5CvqRc2LwE8WQ) will generate the control files for the CodeML runs. <br>
     * Set up bash loop for creating these files based on lists of gene sequence alignment files and the appropriate tree file. <br>
     * Use the [**rooted_tree.nwk**]() file as base for labeled input trees (as above).
         * Need to use tree without "*Root*" node label in CodeML runs, but does not need to be truly unrooted. <br>
  6. Make hypothesis-specific (i.e. subdirectory specific) shell scripts from [**lrt_construction_template.sh**](https://hackmd.io/urgsFStZS4WiIfFI_YbWFA) and [**run_codeml_template.sh**](https://hackmd.io/tn690QwRS5Ca5FzqH1Djrw).
      * These files must be in the hypothesis specific subdirectory to run properly. 
      * For run_codeml_template.sh, must modify NAMESFILE (list of alignment files), and TREE (path to appropriate Newick tree).
      * For lrt_construction_template.sh, must modify job name and codeml_<HYPOTHESIS>_output.txt name.

## CodeML Analysis
  6. Run CodeML by submitting the appropriate run_codeml_<HYPOTHESIS>.sh script from within the hypothesis subdirectory (i.e. selection/long_lived/).
     * Calls ctl_files.py to make appropriate CodeML control files then runs CodeML. 
  8. After array job completion, collate LRT values for each gene in the hypothesis subdirectory by submitting the appropriate lrt_construction_<HYPOTHESIS>.sh script (script must be in hypothesis subdirectory to run properly).
     * Processes and reformats CodeML output.
       * Extracts and concatenates likelihood values from all null and alternative model CodeML runs within the directory's subdirectories, then calls python_test.py to create CSV file of LRT values, then finally calculates the p-value of LRTs. 
     * Final output is codeml_<HYPOTHESIS>_output.txt
<br>

## Other Notes

  * Best set up is likely running CodeML first to find TE defense and DNA repair genes with selection signal in vesper bats, then use the subset with signal for the pairwise comparisons for long- vs short-lived species (if limiting to comparisons within Vespers).

  * NOTE: CodeML can be installed as a python submodule of Biopython.PAML (see doc: (https://biopython.org/wiki/PAML), but unclear how to call outside of python session). <br>


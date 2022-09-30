## Setting up
<br>
  1. Run codeml_input.py to generate the control files for the CodeML runs
  * Set up bash loop for creating these files based on lists of gene sequence alignment files and the appropriate tree file
  * Use the unrooted_tree.nwk file as base (need to use unrooted trees in CodeML runs)
    * To make the unrooted_tree.nwk (or to modify it) in R, use ape
    ```
    library(ape)

    bat_tr=read.tree("03.16.22_tree.tre")
    bat_rooted=(root(bat_tr, 16:17, resolve.root = TRUE))
    is.rooted(bat_rooted)
    plot(bat_rooted)
    tiplabels(cex = .75, bg = "yellow")
    bat2=(unroot(bat_rooted))
    plot(bat2, show.node.label = TRUE)
    bat4=(makeNodeLabel(bat2, c("u"), prefix="#1", nodeList = list(Bat="HL")))
    plot(bat3, show.node.label = TRUE)
    is.rooted(bat2) + nodelabels(frame = "r", bg = "yellow", adj = 0)
    write.tree(bat2, file="unrooted_tree.nwk")
    ```
    
    * Remove the #1 from the ancestral bat branch if need to reroot tree to human/mouse outgroups, then unroot again
    * Move #1 from the ancestral bat branch to the base of Vespertilionidae (can do in tree view to make sure it's on the correct node)
      * This is for selection testing if Vespertilionid bats have a different signal than other bats
    * Move #1 from the ancestral bat branch to the base of each short-lived species (1 per tree, for pairwise comparisons)
      * This is for selection testing to determine if short-lived species have different signal from the long-lived species
 <br>
  2. Run CodeML
  * Best set up is likely running CodeML first to find TE defense and DNA repair genes with selection signal in vesper bats, then use the subset with signal for the pairwise comparisons for long- vs short-lived species (if limiting to comparisons within Vespers)
  * Null model is that all taxa have same selection signal (foreground = background)
  * Alternative model is that different branches have different selection signals (foreground =/= background)
  * CodeML is installed as a python submodule of Biopython.PAML (see doc: https://biopython.org/wiki/PAML)
  <br>
  3. Run likelihood_test.py (which imports lrt_construction.py) to process and reformat CodeML output

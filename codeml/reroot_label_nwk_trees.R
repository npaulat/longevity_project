library(ape)
library(treeio)
library(ggtree)
library(tidytree)
## Also uses functions from Gmisc and xfun ##

#### Set Up ####
setwd("/lustre/scratch/npaulat/selection/trees")

### Get Input Newick Tree from Command Line ###
args <- commandArgs(trailingOnly = TRUE)
## ex: Rscript reroot_label_nwk_trees.R ENST00000374900.FAAH2.cleanLb_hmm_manual.fasta.nwk
## Check for argument
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  # default output file
  stop("Too many arguments supplied, should be one (input file).n", call.=FALSE)
}

tree_file <- args[1]
tree_base <- basename(tree_file)
gene_base <- paste(strsplit(tree_base, "[.]")[[1]][1:2], collapse = ".")


#### Manipulate Input Tree ####

### Reroot input gene tree ###
tree=read.tree(tree_file)
tree$node.label <- NULL
tree_info <- as_tibble(tree)
step1 <- tree_info %>% filter(!grepl('HL|100', label))
out_info <- step1 %>% filter(!is.na(label))
out_nodes <- range(out_info$node)

#bat_rooted=(root(tree, 12:18, resolve.root = TRUE)) # <- for tree with all 30 species ("forcedtop_tree.nwk")
bat_rooted=(root(tree, c(out_nodes[1]:out_nodes[2]), resolve.root = TRUE))
#ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
#ggtree(bat_rooted) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
is.rooted(bat_rooted)
#plot(bat_rooted, show.tip.label = TRUE)
## This appears to be properly equivalent to the topology and branch lengths of "chiroptera_branch.nwk", 
## without node labels
out_dir <- "/lustre/scratch/npaulat/selection/trees/rooted_trees/"
tree_name <- paste(gene_base, ".rooted_tree.nwk", sep = "")
tree_file <- Gmisc::pathJoin(out_dir, tree1_name)
#write.tree(bat1, file="tree_chiroptera_branch.nwk")
write.tree(bat1, file= tree1_file)
###

### Make chiroptera_branch tree (tree 1) ###
bat1=(makeNodeLabel(bat_rooted, method="user", nodeList = list("#1" = c("HL"))))
#plot(bat1, show.node.label = TRUE)
## The next step should be redundant, but just in case
bat1$node.label <- replace(bat1$node.label, c(1:1), "")
#ggtree(bat1) + geom_text2(aes(label=branch.length), hjust=-.1) + geom_tiplab()
out_dir1 <- "/lustre/scratch/npaulat/selection/trees/chiroptera_trees/"
tree1_name <- paste(gene_base, ".chiroptera_branch.nwk", sep = "")
tree1_file <- Gmisc::pathJoin(out_dir, tree1_name)
#write.tree(bat1, file="tree_chiroptera_branch.nwk")
#write.tree(bat1)
write.tree(bat1, file=tree1_file)
###

### Make long_lived_branch tree (tree 2) ###
batNames2 <- data.frame(old.labels = c("HLantPal2", "HLmyoMyo6", "HLtadBra2", "HLdesRot8A", "HLartJam3"),
                                     new.labels = c("HLantPal2 #1", "HLmyoMyo6 #1", "HLtadBra2 #1", "HLdesRot8A #1", "HLartJam3 #1"))
bat2 <- bat1
bat2$node.label <- NULL
bat2$tip.label[bat2$tip.label %in% batNames2$old.labels] <- batNames2$new.labels
bat2$tip.label <- gsub("_", " ", bat2$tip.label)
#plot(bat2, show.node.label = TRUE)
#ggtree(bat2) + geom_text2(aes(label=branch.length), hjust=-.1) + geom_tiplab()
out_dir2 <- "/lustre/scratch/npaulat/selection/trees/long_lived_trees/"
tree2_name <- paste(gene_base, ".chiroptera_branch.nwk", sep = "")
tree2_file <- Gmisc::pathJoin(out_dir, tree2_name)
#write.tree(bat2, file="tree_long_lived_branch.nwk")
#write.tree(bat2)
write.tree(bat2, file=tree2_file)
## Will have to sed change "_" to " " because writing .nwk file converts any spaces to underscores
#xfun::gsub_file("tree_long_lived_branch.nwk", "_#1", " #1", fixed = TRUE)
xfun::gsub_file(tree2_file, "_#1", " #1", fixed = TRUE)
###

### Make short_lived_branch tree (tree 3) ###
batNames3 <- data.frame(old.labels = c("HLmyoNig1", "HLmolMol2", "HLdipEca1", "HLsacBil1"),
                       new.labels = c("HLmyoNig1 #1", "HLmolMol2 #1", "HLdipEca1 #1", "HLsacBil1 #1"))
bat3 <- bat1
bat3$node.label <- NULL
bat3$tip.label[bat2$tip.label %in% batNames3$old.labels] <- batNames3$new.labels
bat3$tip.label <- gsub("_", " ", bat3$tip.label)
#plot(bat3, show.node.label = TRUE)
#ggtree(bat3) + geom_text2(aes(label=branch.length), hjust=-.1) + geom_tiplab()
out_dir3 <- "/lustre/scratch/npaulat/selection/trees/short_lived_trees/"
tree3_name <- paste(gene_base, ".short_lived_branch.nwk", sep = "")
tree3_file <- Gmisc::pathJoin(out_dir, tree3_name)
#write.tree(bat3, file="tree_short_lived_branch.nwk")
#write.tree(bat3)
write.tree(bat3, tree3_file)
## Will have to sed change "_" to " " bc writing .nwk file converts any spaces to underscores
#xfun::gsub_file("tree_short_lived_branch.nwk", "_#1", " #1", fixed = TRUE)
xfun::gsub_file(tree3_file, "_#1", " #1", fixed = TRUE)
###

### For vespertilionidae_branch tree (tree 4) ###
batNames4 <- c("HLantPal2", 'HLpipKuh2', "HLmyoMyo6", "HLmyoNig1")
vesperInfo <- subset(tree_info, (label %in% batNames4))
vesperSp <- vesperInfo$label

bat4=(makeNodeLabel(bat_rooted, method="user", nodeList = list("#1" = vesperSp)))
#plot(bat4, show.node.label = TRUE)
# The next step should be redundant, but just in case
bat4$node.label <- replace(bat4$node.label, c(1:1), "")
#ggtree(bat4) + geom_text2(aes(label=branch.length), hjust=-.1) + geom_tiplab()
out_dir4 <- "/lustre/scratch/npaulat/selection/trees/vesper_trees/"
tree4_name <- paste(gene_base, ".vespertilionidae_branch.nwk", sep = "")
tree4_file <- Gmisc::pathJoin(out_dir, tree4_name)
#write.tree(bat4, file="tree_vespertilionidae_branch.nwk")
write.tree(bat4, file=tree4_file)
###

### For rhinolophidae_branch tree (tree 5) ###
batNames5 <- c("HLrhiFer5", 'HLrhiAff2', "HLrhiSed2", "HLrhiLuc4", "HLrhiTri2")
rhiInfo <- subset(tree_info, (label %in% batNames5))
rhiSp <- rhiInfo$label

out_dir5 <- "/lustre/scratch/npaulat/selection/trees/rhinolophid_trees/"
tree5_name <- paste(gene_base, ".rhinolophidae_branch.nwk", sep = "")
tree5_file <- Gmisc::pathJoin(out_dir, tree5_name)

if (length(rhiSp) > 1) {
  bat5=(makeNodeLabel(bat_rooted, method="user", nodeList = list("#1" = rhiSp)))
  #plot(bat5, show.node.label = TRUE)
  # The next step should be redundant, but just in case
  bat5$node.label <- replace(bat5$node.label, c(1:1), "")
  #ggtree(bat5) + geom_text2(aes(label=branch.length), hjust=-.1) + geom_tiplab()
  #write.tree(bat5, file="tree_rhinolophidae_branch.nwk")
  write.tree(bat5, file=tree5_file)
} else {
  bat5 = bat1
  bat5$node.label <- NULL
  new_label <- paste(rhiSp[1], "#1", sep=" ")
  bat5$tip.label[bat5$tip.label==rhiSp] <- new_label
  #write.tree(bat5, file="tree_rhinolophidae_branch.nwk")
  write.tree(bat5, tree5_file)
  ### Will have to sed change "_" to " " bc writing .nwk file converts any spaces to underscores
  xfun::gsub_file("tree_rhinolophidae_branch.nwk", "_#1", " #1", fixed = TRUE)
}
###

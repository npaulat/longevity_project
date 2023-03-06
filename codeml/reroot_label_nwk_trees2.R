library(ape)
library(treeio)
library(ggtree)
library(tidytree)
library(data.table)

setwd("C:/Users/Nikki/OneDrive - Texas Tech University/longevity_trees/")
setwd("C:/Users/Nikki/Downloads/longevity_trees/")
#setwd("/lustre/scratch/npaulat/selection/trees")

tree_list <- "C:/Users/Nikki/OneDrive - Texas Tech University/longevity_trees/genes_with_min_species.txt"

#tree_file <- "C:/Users/Nikki/OneDrive - Texas Tech University/ENST00000254719.RPA1.cleanLb_hmm_manual.fasta.nwk" # <- all 30 sp

reworkTree <- function(tree_file){

  tree_base <- basename(tree_file)
  #tree_base <- basename("C:/Users/Nikki/OneDrive - Texas Tech University/ENST00000374900.FAAH2.cleanLb_hmm_manual.fasta.nwk")
  gene_base <- paste(strsplit(tree_base, "[.]")[[1]][1:2], collapse = ".")
  
  
  #### Manipulate Input Tree ####
  
  ### Reroot input gene tree ###
  tree=read.tree(tree_file)
  tree$node.label <- NULL
  tree_info <- as_tibble(tree)
  step1 <- tree_info %>% filter(!grepl('HL|100', label))
  out_info <- step1 %>% filter(!is.na(label))
  out_nodes <- range(out_info$node)
  
  #bat_rooted=(root(tree, 12:18, resolve.root = TRUE))
  bat_rooted=(root(tree, c(out_nodes[1]:out_nodes[2]), resolve.root = TRUE))
  ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
  ggtree(bat_rooted) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
  is.rooted(bat_rooted)
  #plot(bat_rooted, show.tip.label = TRUE)
  
  # This appears to be properly equivalent to the topology and branch lengths of "chiroptera_branch.nwk", 
  # without node labels
  #out_dir <- "C:/Users/Nikki/OneDrive - Texas Tech University/longevity_trees/"
  out_dir <- "C:/Users/Nikki/Downloads/longevity_trees/genes_min_sp_labeled/"
  #out_dir <- "/lustre/scratch/npaulat/selection/trees/rooted_trees/"
  tree_name <- paste(gene_base, ".rooted_tree.nwk", sep = "")
  tree_file <- Gmisc::pathJoin(out_dir, tree_name)
  write.tree(bat_rooted, file="rooted_tree.nwk")
  write.tree(bat_rooted, file= tree_file)
  write.tree(bat_rooted)
  ###
  
  ### Make chiroptera_branch tree (tree 1) ###
  #bat1 <- bat_rooted
  #bat1$node.label <- replace(bat1$node.label, c(1:length(bat1$node.label)), "")
  bat1=(makeNodeLabel(bat_rooted, method="user", nodeList = list("#1" = c("HL"))))
  #plot(bat1, show.node.label = TRUE)
  
  # The next step should be redundant, but just in case
  bat1$node.label <- replace(bat1$node.label, c(1:1), "")
  #ggtree(bat1) + geom_text2(aes(label=branch.length), hjust=-.1) + geom_tiplab()
 
   #out_dir1 <- "/lustre/scratch/npaulat/selection/trees/chiroptera_trees/"
  out_dir1 <- out_dir
  tree1_name <- paste(gene_base, ".chiroptera_branch.nwk", sep = "")
  tree1_file <- Gmisc::pathJoin(out_dir1, tree1_name)
  #write.tree(bat1, file="tree_chiroptera_branch.nwk")
  #write.tree(bat1)
  write.tree(bat1, file=tree1_file)
  ###
  
  ### Make long_lived_branch tree (tree 2) ###
  batNames2a <- data.frame(old.labels = c("HLantPal2"),
                          new.labels = c("HLantPal2 #1"))
  batNames2b <- data.frame(old.labels = c("HLmyoMyo6"),
                          new.labels = c("HLmyoMyo6 #1"))
  batNames2c <- data.frame(old.labels = c("HLtadBra2"),
                          new.labels = c("HLtadBra2 #1"))
  batNames2d <- data.frame(old.labels = c("HLdesRot8A"),
                          new.labels = c("HLdesRot8A #1"))
  batNames2e <- data.frame(old.labels = c("HLartJam3"),
                          new.labels = c("HLartJam3 #1"))
  bat2 <- bat_rooted
  bat2$node.label <- NULL
  bat2$tip.label[bat2$tip.label %in% batNames2a$old.labels] <- batNames2a$new.labels
  bat2$tip.label[bat2$tip.label %in% batNames2b$old.labels] <- batNames2b$new.labels
  bat2$tip.label[bat2$tip.label %in% batNames2c$old.labels] <- batNames2c$new.labels
  bat2$tip.label[bat2$tip.label %in% batNames2d$old.labels] <- batNames2d$new.labels
  bat2$tip.label[bat2$tip.label %in% batNames2e$old.labels] <- batNames2e$new.labels
  bat2$tip.label <- gsub("_", " ", bat2$tip.label)
  ##plot(bat2, show.node.label = TRUE)
  #ggtree(bat2) + geom_text2(aes(label=branch.length), hjust=-.1) + geom_tiplab()
  
  #out_dir2 <- "/lustre/scratch/npaulat/selection/trees/long_lived_trees/"
  out_dir2 <- out_dir
  tree2_name <- paste(gene_base, ".long_lived_branch.nwk", sep = "")
  tree2_file <- Gmisc::pathJoin(out_dir2, tree2_name)
  #write.tree(bat2, file="tree_long_lived_branch.nwk")
  #write.tree(bat2)
  write.tree(bat2, file=tree2_file)
  ### Will have to sed change "_" to " " because writing .nwk file converts any spaces to underscores
  #xfun::gsub_file("xx_long_lived_branch.nwk", "_#1", " #1", fixed = TRUE)
  xfun::gsub_file(tree2_file, "_#1", " #1", fixed = TRUE)
  ###
  
  ### Make short_lived_branch tree (tree 3) ###
  batNames3a <- data.frame(old.labels = c("HLmyoNig1"),
                          new.labels = c("HLmyoNig1 #1"))
  batNames3b <- data.frame(old.labels = c("HLmolMol2"),
                          new.labels = c("HLmolMol2 #1"))
  batNames3c <- data.frame(old.labels = c("HLdipEca1"),
                          new.labels = c("HLdipEca1 #1"))
  batNames3d <- data.frame(old.labels = c("HLsacBil1"),
                          new.labels = c("HLsacBil1 #1"))
  bat3 <- bat_rooted
  bat3$node.label <- NULL
  bat3$tip.label[bat3$tip.label %in% batNames3a$old.labels] <- batNames3a$new.labels
  bat3$tip.label[bat3$tip.label %in% batNames3b$old.labels] <- batNames3b$new.labels
  bat3$tip.label[bat3$tip.label %in% batNames3c$old.labels] <- batNames3c$new.labels
  bat3$tip.label[bat3$tip.label %in% batNames3d$old.labels] <- batNames3d$new.labels
  bat3$tip.label <- gsub("_", " ", bat3$tip.label)
  #plot(bat3, show.node.label = TRUE)
  #ggtree(bat3) + geom_text2(aes(label=branch.length), hjust=-.1) + geom_tiplab()
  
  #out_dir3 <- "/lustre/scratch/npaulat/selection/trees/short_lived_trees/"
  out_dir3 <- out_dir
  tree3_name <- paste(gene_base, ".short_lived_branch.nwk", sep = "")
  tree3_file <- Gmisc::pathJoin(out_dir3, tree3_name)
  #write.tree(bat3, file="tree_short_lived_branch.nwk")
  #write.tree(bat3)
  write.tree(bat3, tree3_file)
  ### Will have to sed change "_" to " " bc writing .nwk file converts any spaces to underscores
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
  
  #out_dir4 <- "/lustre/scratch/npaulat/selection/trees/vesper_trees/"
  out_dir4 <- out_dir
  tree4_name <- paste(gene_base, ".vespertilionidae_branch.nwk", sep = "")
  tree4_file <- Gmisc::pathJoin(out_dir4, tree4_name)
  #write.tree(bat4, file="tree_vespertilionidae_branch.nwk")
  write.tree(bat4, file=tree4_file)
  ###
  
  ### For rhinolophidae_branch tree (tree 5) ###
  batNames5 <- c("HLrhiFer5", 'HLrhiAff2', "HLrhiSed2", "HLrhiLuc4", "HLrhiTri2")
  rhiInfo <- subset(tree_info, (label %in% batNames5))
  rhiSp <- rhiInfo$label
  
  #out_dir5 <- "/lustre/scratch/npaulat/selection/trees/rhinolophid_trees/"
  out_dir5 <- out_dir
  tree5_name <- paste(gene_base, ".rhinolophidae_branch.nwk", sep = "")
  tree5_file <- Gmisc::pathJoin(out_dir5, tree5_name)
  
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
    #xfun::gsub_file("tree_rhinolophidae_branch.nwk", "_#1", " #1", fixed = TRUE)
    xfun::gsub_file(tree5_file, "_#1", " #1", fixed = TRUE)
  }
}
#|>
###

gene_list <- read.table(tree_list, header=F)[[1]]
for (k in 1:length(gene_list)) {
  gene_file <- paste("C:/Users/Nikki/Downloads/longevity_trees/genes_min_sp_unlabeled/", gene_list[k], ".nwk", sep="")
  reworkTree(gene_file)
}

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if (!require("pacman")) install.packages("pacman")

pacman::p_load(ggtree, ggimage)

newick.tree <- read.tree(args[1])

newick.tree$tip.label <- gsub("_", " ", newick.tree$tip.label)

#phylopic images
tips <- newick.tree$tip.label
tips[6] <- "Eptesicus"
tips[8] <- "Marmota monax"
tips[10] <- "Mus musculus"
tipsimg <- ggimage::phylopic_uid(tips)
tipsimg$name <- tips
#save time and also edit if necessary
#write.csv(tipsimg, "tipsimg.csv")
#tipsimg <- read.csv("tipsimg.csv")

ggtree(newick.tree) + 
  geom_label2(aes(subset=!isTip, label=label), col="red4", size=2, fill="white", label.size=0.2) + 
  geom_tiplab(image=tipsimg$uid, geom="phylopic", offset = .075) +
  xlim(0,.4) + 
  geom_treescale(x=0, y=11) +
  geom_tiplab(fontface="italic")

ggsave(args[2], width=12, height=6, dpi=600)

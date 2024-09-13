#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if (!require("pacman")) install.packages("pacman")

pacman::p_load(ggtree)

newick.tree <- read.tree(args[1])

newick.tree$tip.label <- gsub("_", " ", newick.tree$tip.label)

ggtree(newick.tree) + 
  geom_label2(aes(subset=!isTip, label=label), col="red4", size=2, fill="white", label.size=0.2) + 
  xlim(0,.4) + 
  geom_treescale(x=0, y=11) + 
  geom_tiplab(aes(image=paste0(args[2], label, '.jpeg')), 
            geom="image", offset=.1, align=TRUE, size=.1, linesize = 1) +
  geom_tiplab(geom="label", fontface="italic", fill="white", alpha=1)

ggsave(args[3], width=12, height=6, dpi=600)

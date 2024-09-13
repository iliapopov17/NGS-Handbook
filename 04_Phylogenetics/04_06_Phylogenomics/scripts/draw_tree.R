#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(ggtree)

tr <- read.tree(args[1])

ggtree(tr) + geom_tiplab() + 
  xlim(0,1) + geom_label2(aes(subset=!isTip, label=label), col="red4", size=2, fill="white", label.size=0.2) + 
  geom_treescale()

ggsave(args[2], width = 8, height=6)
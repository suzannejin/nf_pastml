#!/usr/bin/env Rscript

library('ape')

args = commandArgs(trailingOnly=TRUE)
treename = args[1]
outname = args[2]

# read tree
tree = read.tree(treename)

# remove branch length
tree$edge.length = NULL

# write modified tree
write.tree(tree, file=outname)
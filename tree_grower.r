#!/usr/bin/env R
# This script generates a phylogenetic tree plot with sequence alignment using ggtree and tidytree.


require(treeio)
require(ggtree)
require(tidytree)
require(data.table)
require(ggplot2)

name = commandArgs(trailingOnly=TRUE)
if (length(name)==0 | length(name)>1) {
  stop("One argument must be supplied (input file).n", call.=FALSE)
}

treefile = "in-tree.txt"
fastafile = "in-fasta.txt"
colors = c("a" = "chartreuse2", "t" = "coral2", "c" = "deepskyblue1", "g" = "darkgoldenrod2", "-" = "white")
treeplot = read.newick(treefile) %>% ggtree()
duoplot = msaplot(treeplot, fastafile, bg_line = FALSE) +
        scale_fill_manual(values = colors,
                          name = "Nucleotide") +
        scale_y_continuous(expand=c(0, 0.3)) +
        labs(title = name, subtitle = "IAPLTR1 clustering and sequence alignment")
ggsave("outfile.png", duoplot, device = "png")


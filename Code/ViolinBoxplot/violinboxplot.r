#!/usr/bin/env Rscript

library(ggplot2)

coprarna_PCC <- read.table('OrdinalRanked_Final_CopraRNA_PCC.txt', header = TRUE, sep = "\t")
snrarftarget_PCC <- read.table('OrdinalRanked_Final_sRNARFTarget_PCC.txt', header = TRUE, sep = "\t")
intarna_PCC <- read.table('OrdinalRanked_Final_IntaRNA_PCC.txt', header = TRUE, sep = "\t")

tmp <- rbind(coprarna_PCC, snrarftarget_PCC, intarna_PCC)
Cor <- tmp$OrdinalRank
Program <- tmp$Program

ggplot(tmp, aes(y=Cor, x=Program, fill = Program)) +  geom_violin() + ylim(1,6344) + ylab("Ranking of True Positives") + geom_boxplot(width=0.1) +theme(legend.position = "top")


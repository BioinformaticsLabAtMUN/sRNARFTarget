#!/usr/bin/env Rscript

library(randomForest)
library(PRROC)

res102 <- read.table('./Data/PlotsData/PRROC/Ecoli/sRNARFTarget_ecoli_class1.txt', header = FALSE, sep = "\t")
resrem <- read.table('./Data/PlotsData/PRROC/Ecoli/sRNARFTarget_ecoli_class0.txt', header = FALSE, sep = "\t")


Pred102 <- res102[,3]
Predrem <- resrem[,3]

label102 <- res102[,4]
labelrem <- resrem[,4]

x<-c(Pred102,Predrem);
lab<-c(label102,labelrem)

roc <- roc.curve(scores.class0 = x, weights.class0 = lab, curve = TRUE);
plot(roc);
print(roc);

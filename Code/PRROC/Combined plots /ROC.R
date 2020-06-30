#!/usr/bin/env Rscript

library(randomForest)
library(PRROC)

res22 <- read.table('CopraRNA_ecoli_class1.txt', header = FALSE, sep = "\t")
resrem <- read.table('CopraRNA_ecoli_class0.txt', header = FALSE, sep = "\t")

Pred22 <- res22[,3]
Predrem <- resrem[,3]

label22 <- res22[,4]
labelrem <- resrem[,4]

x<-c(Pred22,Predrem);
lab<-c(label22,labelrem)


res22c <- read.table('sRNARFTarget_ecoli_class1.txt', header = FALSE, sep = "\t")
resremc <- read.table('sRNARFTarget_ecoli_class0.txt', header = FALSE, sep = "\t")

Pred22c <- res22c[,3]
Predremc <- resremc[,3]

label22c <- res22c[,4]
labelremc <- resremc[,4]

xc<-c(Pred22c,Predremc);
labc<-c(label22c,labelremc)


res22cc <- read.table('IntaRNA_ecoli_class1.txt', header = FALSE, sep = "\t")
resremcc <- read.table('IntaRNA_ecoli_class0.txt', header = FALSE, sep = "\t")

Pred22cc <- res22cc[,3]
Predremcc <- resremcc[,3]

label22cc <- res22cc[,4]
labelremcc <- resremcc[,4]

xcc<-c(Pred22cc,Predremcc);
labcc<-c(label22cc,labelremcc)


roc1 <- roc.curve(scores.class0 = x, weights.class0 = lab, curve = TRUE, rand.compute = T);

roc2 <- roc.curve(scores.class0 = xc, weights.class0 = labc, curve = TRUE, rand.compute = T);

roc3 <- roc.curve(scores.class0 = xcc, weights.class0 = labcc, curve = TRUE, rand.compute = T);


# plot PR curve in red, without legend
plot( roc1, color = 2, auc.main=TRUE, rand.plot = TRUE);
# add second PR curve in green
plot( roc2, color = 3, add = TRUE);
# add third plot 
plot( roc3, color = 4, add = TRUE);

legend("bottomright", legend=c("CopraRNA", "sRNARFTarget", "IntaRNA"), col=c("red", "green", "blue"), lty = c(1, 1), lwd = c(2,2),cex=0.8, box.lty=0)

print(roc1);
print(roc2);
print(roc3);





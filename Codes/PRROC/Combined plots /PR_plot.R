#!/usr/bin/env Rscript

library(randomForest)
library(PRROC)

res102 <- read.table('./Data/PlotsData/PRROC/Ecoli/CopraRNA_ecoli_class1.txt', header = FALSE, sep = "\t")
resrem <- read.table('./Data/PlotsData/PRROC/Ecoli/CopraRNA_ecoli_class0.txt', header = FALSE, sep = "\t")


Pred102 <- res102[,3]
Predrem <- resrem[,3]

label102 <- res102[,4]
labelrem <- resrem[,4]

x<-c(Pred102,Predrem);
lab<-c(label102,labelrem)


res22c <- read.table('./Data/PlotsData/PRROC/Ecoli/sRNARFTarget_ecoli_class1.txt', header = FALSE, sep = "\t")
resremc <- read.table('./Data/PlotsData/PRROC/Ecoli/sRNARFTarget_ecoli_class0.txt', header = FALSE, sep = "\t")

Pred22c <- res22c[,3]
Predremc <- resremc[,3]

label22c <- res22c[,4]
labelremc <- resremc[,4]

xc<-c(Pred22c,Predremc);
labc<-c(label22c,labelremc)


res22cc <- read.table('./Data/PlotsData/PRROC/Ecoli/IntaRNA_ecoli_class1.txt', header = FALSE, sep = "\t")
resremcc <- read.table('./Data/PlotsData/PRROC/Ecoli/IntaRNA_ecoli_class0.txt', header = FALSE, sep = "\t")

Pred22cc <- res22cc[,3]
Predremcc <- resremcc[,3]

label22cc <- res22cc[,4]
labelremcc <- resremcc[,4]

xcc<-c(Pred22cc,Predremcc);
labcc<-c(label22cc,labelremcc)


pr1 <- pr.curve(scores.class0 = x, weights.class0 = lab, curve = TRUE, rand.compute = T);

pr2 <- pr.curve(scores.class0 = xc, weights.class0 = labc, curve = TRUE, rand.compute = T);

pr3 <- pr.curve(scores.class0 = xcc, weights.class0 = labcc, curve = TRUE, rand.compute = T);

# plot PR curve in red, without legend 2 for red 3 for green
#plot( pr1, color = 2, auc.main="FALSE",main = "Ecoli",rand.plot = TRUE );
plot( pr1, color = 2, auc.main="FALSE",main = expression(italic("Escherichia coli")),rand.plot = TRUE );

# add second PR curve in green
plot( pr2, color = 3, add = TRUE );

plot( pr3, color = 4, add = TRUE );

legend("topright", legend=c("CopraRNA", "sRNARFTarget", "IntaRNA"), col=c("red", "green", "blue"), lty = c(1, 1), lwd = c(2,2),cex=0.8, box.lty=0)

print(pr1);
print(pr2);
print(pr3);






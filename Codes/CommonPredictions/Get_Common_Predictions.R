#!/usr/bin/env Rscript

#Read the files
RF <- read.table("./Data/ProgramResults/OriginalPredictions/Ecoli/Final_sRNARFTarget_Ecoli.txt", header = T, sep= "\t",stringsAsFactors = F)
IntaRNA <- read.table("./Data/ProgramResults/OriginalPredictions/Ecoli/Final_IntaRNA_Ecoli.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(RF)
dim(IntaRNA)

# Get the rows in common
inCommon <- merge(IntaRNA, RF, by = c("sRNA", "mRNA"), sort = T)
dim(inCommon)
head(inCommon)

#Create two files with the rows in common

IntaRNA_inCommon <- inCommon[,c("sRNA", "mRNA", "IntaRNA_Score", "Class.x")]
colnames(IntaRNA_inCommon) <- c("sRNA", "mRNA", "IntaRNA_Score", "Class")

write.csv(IntaRNA_inCommon, file = "IntaRNA_Ecoli_inCommon_IR.csv")

RF_incommon <- inCommon[,c("sRNA", "mRNA", "RF_Score","Class.y")]
colnames(RF_incommon) <- c("sRNA", "mRNA", "RF_Score", "Class")

write.csv(RF_incommon, file = "RF_Ecoli_incommon_IR.csv")

dim(IntaRNA_inCommon)
dim(RF_incommon)

#!/usr/bin/env Rscript

r <- read.table('testdf.txt', header = T, sep = "\t")

selectMostSignificant <- function(r){
     x <- by(r, cbind(r['sRNA'], r['mRNA']), function(m){ m[which(m[,"pvalue"] ==
min(m[,"pvalue"]))[1], , drop = FALSE  ]  }, simplify = FALSE)
     x <- do.call("rbind", x)
     colnames(x) <- c("sRNA", "mRNA", "pvalue")
     x
}

result <- selectMostSignificant(r)


write.table(result, file = multocidasignificantonesresult.txt”, sep = "\t")


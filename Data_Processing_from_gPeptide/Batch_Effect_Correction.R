if(!require(sva))
  BiocManager::install("sva")
if(!require(bladderbatch))
  BiocManager::install("bladderbatch")

library(sva)
library(bladderbatch)

em <- read.csv('Expression_Matrix_Imputed.csv', header = TRUE, row.names = 1)
batch <- c(3, 1, 3, 2, 2, 1, 1, 3, 1, 2, 3, 3, 3, 3, 2, 3, 2, 3)
combat_em <- ComBat(dat=em, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

write.csv(combat_em, 'Expression_Matrix_BEC.csv')

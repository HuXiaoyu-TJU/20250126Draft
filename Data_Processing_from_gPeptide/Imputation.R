library(missForest)

set.seed(0)
em <- read.csv('Expression_Matrix.csv', header = TRUE, row.names = 1)
em_rf <- missForest(em,ntree = 500)
em_rf <- em_rf$ximp
write.csv(em_rf, 'Expression_Matrix_Imputed.csv')

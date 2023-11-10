
######alpha_diversity
library(vegan)
library(picante)      
setwd("")
alpha_diversity <- function(x, tree = NULL) {
  Chao1 <- estimateR(x)[2, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Shannon <- sprintf("%0.4f", Shannon)
  result <- data.frame(Chao1, Shannon)
}
otu <- read.delim('otu.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)
alpha <- alpha_diversity (otu)
write.csv(alpha, 'alpha_diversity.csv', quote = FALSE)

library(ConsensusClusterPlus)        
data=read.table("genus.txt" , header=T, sep="\t", check.names=F, row.names=1)
data <- data[rowSums(data == 0) < 0.4* ncol(data),]
dim(data)
data <- log2(data+1)
data = t(scale(t(data)))
data=as.matrix(data)
dim(data)
maxK=4
results=ConsensusClusterPlus(data,
              maxK=maxK,
              pItem=0.8,
              pFeature=1,
              reps =1000, 
              clusterAlg="pam",
              distance="spearman",
              seed=123456,
              plot="png")
cluster=cbind(results[[2]][["consensusClass"]],results[[3]][["consensusClass"]],results[[4]][["consensusClass"]])
write.table(cluster,file="cluster.txt",sep="\t",quote=F,col.names=F)

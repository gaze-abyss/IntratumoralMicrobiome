work_dir <- ''
########################
library(igraph)
library(ggsci)
library(scales)
library(stringi)
library(stringr)

get_col <- function(one){
  one_list <- names(table(one))
  all_col <- pal_nejm()(10)[4:8][1:length(one_list)]
  names(all_col) <- one_list
  res <- c()
  for(i in one){
    res <- c(res, all_col[i])
  }
  return(res)
}
work_dir <- ''
G <- list()
set.seed(3)

exp1 <- read.delim(paste0(work_dir,'group1.txt'), row.names = 1, sep = '\t', check.names = FALSE)
pvals1 <- read.delim(paste0(work_dir,'boot1/pvals.two_sided.txt'), row.names = 1, sep = '\t', check.names = FALSE)
corm1 <- read.delim(paste0(work_dir,'cor_sp_t1.txt'), row.names = 1, sep = '\t', check.names = FALSE)

tmp = rownames(exp1)
tmp = str_split(tmp,"g__",simplify = T,n = 2)
rownames(exp1) = tmp[,2]
tmp = colnames(pvals1)
tmp = str_split(tmp,"g__",simplify = T,n = 2)
colnames(pvals1) = tmp[,2]
tmp = colnames(corm1)
tmp = str_split(tmp,"g__",simplify = T,n = 2)
colnames(corm1) = tmp[,2]

cnet1 <- graph_from_adjacency_matrix(as.matrix(corm1),weight=T)

e_cnet1 <- as_data_frame(cnet1,'edge')
w_list1 <- e_cnet1$weight
names(w_list1) <- paste0(e_cnet1[,1],'-',e_cnet1[,2])

corm1[corm1 <= 0] <- -1
corm1[corm1 > 0] <- 1
corm1[corm1 == -1] <- 0
corm1 <- 1

pvals1[pvals1 >= 0.05] <- -1
pvals1[pvals1 < 0.05 & pvals1 >= 0] <- 1
pvals1[pvals1 == -1] <- 0
adj1 <- pvals1 * corm1
diag(adj1) <- 0

adj <- adj1

net <- graph_from_adjacency_matrix(as.matrix(adj), mode = 'undirected')
l <- layout_with_graphopt(net,charge=0.1,niter = 5000)
net <- graph_from_adjacency_matrix(as.matrix(adj1), mode = 'undirected')
clu <- components(net)
conn <- groups(clu)
big <- -Inf
pick <- 1
for(i in 1:length(conn)){
  if(length(conn[[i]]) > big){
    big <- length(conn[[i]])
    pick <- i
  }
}
V(net)$vertex.color <- '#E18727FF'
edge_table <- as_data_frame(net, what = c("edges"))
edge_table$cor <- unlist(lapply(as.list(w_list1[paste0(edge_table[,1],'-',edge_table[,2])]),function(x) {if(x>0){'+'}else{'-'}} ))
write.table(edge_table,file=paste0(work_dir,'edge_table.1.txt'),quote=F,sep='\t')
E(net)$edge.color <- unlist(lapply(as.list(w_list1[paste0(edge_table[,1],'-',edge_table[,2])]),function(x) {if(x>0){'#B33B2A'}else{'#106CA9'}} ))
abundance_score_1 <- as.numeric(rowMeans(exp1))
abundance_score_mat_1 <- data.frame(
  value = abundance_score_1,
  name = 'abundance_score_1'
)
hub_score_1 <- hub_score(net,scale=F)$vector
score_mat_1 <- data.frame(
  value = hub_score_1,
  name = 'hub_score_1'
)
G[[1]] <- net

abundance_all_score_mat <- rbind(abundance_score_mat_1)
abundance_all_score_mat[,1] <- (abundance_all_score_mat[,1]-min(abundance_all_score_mat[,1]))/(max(abundance_all_score_mat[,1])-min(abundance_all_score_mat[,1]))
abundance_all_score_mat[,1][which(abundance_all_score_mat[,1] < 0.05)] <- 0.05
abundance_score_1 <- abundance_all_score_mat[abundance_all_score_mat[,2] == 'abundance_score_1',][,1]

all_score_mat <- rbind(score_mat_1,score_mat_2,score_mat_3)
all_score_mat[,1] <- (all_score_mat[,1]-min(all_score_mat[,1]))/(max(all_score_mat[,1])-min(all_score_mat[,1]))
write.table(all_score_mat,file=paste0(work_dir,'all_score_mat.txt'),quote=F,sep='\t')

ran <- quantile(as.numeric(abundance_score_1), probs = c(0.01,0.99))
all_score_mat[,1][which(all_score_mat[,1] < 0.2)] <- 0.2

hub_score_1 <- all_score_mat[all_score_mat[,2] == 'hub_score_1',][,1]

V(G[[1]])$vertex.color <- circlize::colorRamp2(c(ran[1],ran[2]),c("pink","#B33B2A"))(abundance_score_1*2)

pdf(paste0(work_dir,'network.all.1.pdf'),width = 20,height = 10)
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(G[[1]],layout=l, 
     edge.color = adjustcolor(E(G[[1]])$edge.color,alpha.f = .5),
     edge.width = 3,
     vertex.color = V(G[[1]])$vertex.color,
     vertex.frame.color = V(G[[1]])$vertex.color,
     vertex.label.color = 'black', 
     vertex.label.cex = 1,
     vertex.size = hub_score_1*15)
text(x = -1, y = 1.1, 'A', cex=2.5)
legend("bottomleft", legend = c("     0.9","     0.6","     0.3","     0"),title.adj=0.5,cex=1,bty='o',col='#E18727FF',y.intersp = 3,pch = 19, pt.cex = (1*8)*seq(1,0.2,-0.2),pt.lwd = 1,title = "Hub centrality score", trace=F)
legend("bottomright", legend = c("Positive correlation","Negative correlation"),lty = 1,lwd=3.5,pt.cex=3.5,y.intersp = 2,col=c('#B33B2A','#106CA9'))
dev.off()
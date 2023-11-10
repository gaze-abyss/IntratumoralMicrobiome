library(vegan)
library(ggplot2)
library(pairwiseAdonis)
setwd("")
otu_raw <- read.table(file="otu.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
otu <- t(otu_raw)
otu.distance <- vegdist(otu)
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
group = colnames(otu_raw)
group = data.frame(group)
colnames(group) = "samples"
group$group = "T"
tmp = grep("N",group$samples)
group$group[tmp] = "N"
df <- merge(pc12,group,by="samples")
color=c("#106CA9",	"#B33B2A","#21804C")
p1<-ggplot(data=df,aes(x=V1,y=V2,
                       color=group,shape=group))+
  theme_bw()+
  geom_point(size=1.8)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"))+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#1597A5","#FFC24B"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank(),legend.position = "top")


p <- p1 + stat_ellipse(data=df,geom = "polygon",level=0.9,linetype = 2,size=0.5,aes(fill=group),alpha=0.2,show.legend = T)
dune.div <- adonis2(otu ~ group, data = group, permutations = 999, method="bray")
dune.div

# install.packages('RobustRankAggreg')
# library("cosmosR")
# library("TopKLists")
library("RobustRankAggreg")
library(dplyr)
setwd("E:\\2021-10-排序算法\\our data\\RankAggregate")
rm(list =ls())
data<-read.table("2021-12-Rank-all.csv",header = T,sep = ",")

stars<-data[rowSums(is.na(data))<4,]
star_glist <- list(as.character(unlist(arrange(na.omit(stars[,c(1,2)]),Rank.mRNA)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,3)]),Rank.Pr)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,4)]),Rank.Phos)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,5)]),Rank.Acety)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,6)]),Rank.path)[1])))

#进行rra
star_rra <- aggregateRanks(star_glist)
#补上rra之后的排名
star_rra <- cbind(star_rra, "rra_rank" = seq(1,nrow(star_rra)))
#可以看到rra针对各个榜单的排名重新给每个明星计算了分数。按照分数就可以获得新的整合后的排名
#star_rra <- left_join(star_rra, stars, by = "Symbol"
m<-merge(star_rra,stars,by.x = "Name",by.y = "Symbol")
write.table(m,file="2021-12-17-Rank-atleast2-final.csv",sep=",")
#################
rm(list =ls())
data<-read.table("2021-12-Rank-onlyp&a.csv",header = T,sep = ",")

stars<-data[rowSums(is.na(data))<4,]
star_glist <- list(as.character(unlist(arrange(na.omit(stars[,c(1,2)]),Rank.mRNA)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,3)]),Rank.Pr)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,4)]),Rank.Phos)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,5)]),Rank.Acety)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,6)]),Rank.path)[1])))

#进行rra
star_rra <- aggregateRanks(star_glist)
#补上rra之后的排名
star_rra <- cbind(star_rra, "rra_rank" = seq(1,nrow(star_rra)))
#可以看到rra针对各个榜单的排名重新给每个明星计算了分数。按照分数就可以获得新的整合后的排名
#star_rra <- left_join(star_rra, stars, by = "Symbol"
m<-merge(star_rra,stars,by.x = "Name",by.y = "Symbol")
write.table(m,file="2021-12-17-Rank-onlyp&a-atleast2.csv",sep=",")

#################
rm(list =ls())
data<-read.table("2021-12-Rank-onlyp&a-nomRNA.csv",header = T,sep = ",")

stars<-data[rowSums(is.na(data))<3,]
star_glist <- list(#as.character(unlist(arrange(na.omit(stars[,c(1,2)]),Rank.mRNA)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,2)]),Rank.Pr)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,3)]),Rank.Phos)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,4)]),Rank.Acety)[1])),
                   as.character(unlist(arrange(na.omit(stars[,c(1,5)]),Rank.path)[1])))

#进行rra
star_rra <- aggregateRanks(star_glist)
#补上rra之后的排名
star_rra <- cbind(star_rra, "rra_rank" = seq(1,nrow(star_rra)))
#可以看到rra针对各个榜单的排名重新给每个明星计算了分数。按照分数就可以获得新的整合后的排名
#star_rra <- left_join(star_rra, stars, by = "Symbol"
m<-merge(star_rra,stars,by.x = "Name",by.y = "Symbol")
write.table(m,file="2021-12-17-Rank-onlyp&a-nomRNA-atleast2.csv",sep=",")

####################################################################
rm(list=ls())

library(pheatmap) 


rm(list=ls())
a<-read.csv("2021-12-17-Rank-onlyp&a-nomRNA-atleast2.csv",header=T)
rownames(a)<-a[,1]
a<-a[,-(1:3)]
a<-a/3000
a[is.na(a)]<-1
#a<-a[,-(1:2)]

a<-as.matrix(a)
annotation_col = data.frame(
  Sample = factor(rep(c("mRNA", "Pr","Prp","Pra","Path"), c(1,1,1,1,1))))####加样本颜色标签
rownames(annotation_col)=colnames(a)
annotation_col = data.frame(
  Sample = factor(rep(c("Pr","Prp","Pra","Path"), c(1,1,1,1))))####加样本颜色标签
rownames(annotation_col)=colnames(a)
a1<-a[1:50,]
ann_colors = list(Sample = c(mRNA = "SkyBlue", Pr = "Green3", Prp = "DarkGoldenrod4", Pra = "Orchid", Path = "RoyalBlue1"))
#plot(1:5,1:5)
ann_colors = list(Sample = c( Pr = "Green3", Prp = "DarkGoldenrod4", Pra = "Orchid", Path = "RoyalBlue1"))

pheatmap(a1,scale = "none", cluster_rows = F,cluster_cols = F,
         clustering_distance_col = "binary", gaps_col = c(1:4),
         annotation_col = annotation_col,annotation_colors = ann_colors,
         #color = cols,
         colorRampPalette(c("OrangeRed","DarkGoldenrod2","RosyBrown3","white"))(100),#color = cols,
         #breaks = seq(0,1,0.01),
         #legend_breaks=seq(0,2000,500),
         fontsize=8, fontsize_row=8,cellwidth = 8, cellheight =8) 

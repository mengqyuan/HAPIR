####Figure 5A
library("pheatmap")
all=all[c(1,2,4,6,3,5),]
bk=c(seq(0.4,1,by=0.1))
pheatmap(as.matrix(all),cluster_rows =F,cluster_cols =F,fontsize=6)
use=colMeans(all)
barplot(height=use,space=1.3)


####Figure 5B
library(ggplot2)
library(ggradar)
colors <-c("#94B3B9", "#EDD69A", "#788816","#C490AA","#F07050","#97A0CB")
p=ggradar(plot.data=all,group.point.size = 3,group.line.width = 1,legend.position="right",
values.radar =c(0,0.5,1))+
scale_color_manual(values = colors)

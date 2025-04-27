#####Figure 4A
out=c()
i=c("Melanoma-GSE100797","Melanoma-GSE115821","NSCLC_GSE135222","Melanoma-PRJEB23709","STAD-PRJEB25780")

auclist=list()
auc=c()
for(j in i){
name=paste(,j,"_200_logit_aucell_hallmark_score.txt",sep="")
pred=read.table(name,sep="\t",header=T)
library(pROC)
result=pred[,1:2]
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc=c(auc,round(auc(rocobj),4))
roc_smooth <- smooth(rocobj) 
auclist[[j]]=rocobj
#auclist[[j]]=roc_smooth
}

co=c("#EDD69A", "#788816","#F07050","#C490AA","#97A0CB")

plot(auclist[[1]],col=co[1],lwd = 3)
plot(auclist[[2]],col=co[2],add=TRUE,lwd = 3)
plot(auclist[[3]],col=co[3],add=TRUE,lwd = 3)
plot(auclist[[4]],col=co[4],add=TRUE,lwd = 3)
plot(auclist[[5]],col=co[5],add=TRUE,lwd = 3)
legend("bottomright", legend = paste(i,auc, sep = ": "),
col = co, lwd = 2, bty = "n",cex = 0.8)


#######Figure 4D
ggplot(re1,aes(x=sig,y=AUC))+geom_boxplot()+geom_point(aes(color = dataset,fill=dataset),
    position = position_jitterdodge(dodge.width = .6, seed = 1),shape=21,
     stroke = 1, size = 3, show.legend = T,
  )+scale_fill_manual(values = c("#94B3B9", "#EDD69A", "#788816","#F07050","#C490AA","#97A0CB"))+scale_color_manual(values = c("#0A7580","#F7834C","#1B2C3E","#DB3C2B","#9B5276","#6F75B0"))+
  theme_classic2()+theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle=60, hjust = 1,vjust=1))+theme(axis.text.x = element_text(color = col))
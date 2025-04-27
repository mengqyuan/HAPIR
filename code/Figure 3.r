######Figure 3A auc
i=c("Melanoma-GSE91061")
score=getexp(j)
PD1=score$exp[which(rownames(score$exp)=="PDCD1"),]
PDL1=score$exp[which(rownames(score$exp)=="CD274"),]
pred=as.data.frame(cbind(as.numeric(PD1),1-as.numeric(score$label)))
rocobj1 <- roc(pred[,2],pred[,1],direction="<")
pred=as.data.frame(cbind(as.numeric(PDL1),1-as.numeric(score$label)))
rocobj2 <- roc(pred[,2],pred[,1],direction="<")

setwd("/mnt/yuanmengqin/Immune_therapy_predict/score_auc/test/10bei")
pred=read.table("Melanoma-GSE91061_200_logit_aucell_hallmark_pred.txt",sep="\t",header=F)
library(pROC)
result=pred[,2:3]
result[,1] <- as.numeric(result[,1])
rocobj3 <- roc(result[,2],result[,1],direction="<")
auc=c()
auc=c(auc,round(auc(rocobj1),4))
auc=c(auc,round(auc(rocobj2),4))
auc=c(auc,round(auc(rocobj3),4))

co=c("#EDD69A", "#788816","#F07050","#C490AA","#97A0CB")
i=c("PD1","PDL1","hallmark_logit")

plot(rocobj1,col=co[1],lwd = 3)
plot(rocobj2,col=co[2],add=TRUE,lwd = 3)
plot(rocobj3,col=co[3],add=TRUE,lwd = 3)
legend("bottomright", legend = paste(i,auc, sep = ": "),
col = co, lwd = 2, bty = "n",cex = 0.8)


#######Figure 3B-D
library("survival")
library("survminer")
merge_data<-data.frame(time=as.numeric(info$survival[,2]),event=as.numeric(os_status),riskScore=as.numeric(pred[rownames(info$survival),2]))
print(head(merge_data))
res_cut <- surv_cutpoint(merge_data, time = "time", event = "event",variables = c("riskScore"))
print(res_cut)
res_cat <- surv_categorize(res_cut)
head(res_cat)
fit <- surv_fit(Surv(time, event) ~riskScore, data = res_cat)
gg1<-ggsurvplot(fit,pval = TRUE,pval.coord = c(100, 0.6),palette =c("#F7834C","#B2C55C"),
          ,ylim = c(0.1, 1), xlab = "Time(Days)",
          risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           linetype = c(1,1), # Change line type by groups
            ggtheme = theme_classic(),
  font.main = c(16),       
  font.x = c(16),              
  font.y = c(16),               
  font.tickslab = c(16),     
  font.legend = c(16),      
  font.caption = c(16),      
  font.risk.table = c(16) 
           )

####Figure 3E-F
library("ggplot2")
p9<-ggplot(final,aes(x=sig,ymin=0,ymax=AUC,group=dataset)) +
 geom_linerange(position = position_dodge2(.5),linetype="dotdash",linewidth=0.5)+
geom_point(aes(x=sig, y=AUC,color=dataset),size=3,position = position_dodge(.5), fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1) +
scale_color_manual(values=c("#2F9894","#D24625"))+
theme_light()+theme(panel.grid.major.x=element_blank(),
panel.border=element_blank(),
axis.ticks.x=element_blank())+
xlab("")+
ylab("AUC")
#scale_y_break(c(2.5,5),space=0.1,scales=0.2)
#coord_flip()
#p9+coord_flip()
p=p9+ theme_classic() +theme(axis.text.x = element_text(angle=45, hjust = 1,vjust=1))+theme(axis.text.x = element_text(color = col))

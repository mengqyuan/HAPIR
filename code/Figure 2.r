#####Figure 2A
dataset="Melanoma-GSE91061"
infile=paste(dataset,"_diff.txt",sep="")
data=read.table(infile,sep="\t")
pvalue=0.05
library(dplyr)
top_up5=top_n(data[which(data[,4]< pvalue),],n=num,wt=logFC)
top_d5=top_n(data[which(data[,4]< pvalue),],n=-num,wt=logFC)

data$sig[!is.element(rownames(data),c(rownames(top_up5),rownames(top_d5)))] <- "NotSig"
data$sig[is.element(rownames(data),rownames(top_up5))] <- "Up"
data$sig[is.element(rownames(data),rownames(top_d5))] <- "Down"
data=data[order(data[,1],decreasing=T),]
data$Marker =0

data1=read.table("Melanoma-GSE91061_200_hallmark_sig.txt",sep="\t")
gene=unique(unlist(strsplit(data1[, 2], "/")))

setwd("/mnt/yuanmengqin/Immune_therapy_predict/score_auc")
data$Marker[is.element(rownames(data),gene)]=1
data$label=ifelse(data$Marker == 1, as.character(rownames(data)), '')
color=c("#234B44","#A6A873","#E28A64")
library(ggrepel)
library(ggplot2)
p=ggplot(data,aes(data$logFC,-1*log10(data$P.Value))) +    
  geom_point(aes(color = sig)) +                           
  labs(title="volcanoplot",                                
       x="log2(FC)", 
       y="-log10(p_val)") + 
  scale_color_manual(values = color) + 
  geom_hline(yintercept=-log10(0.05),linetype=2)+        
  #geom_vline(xintercept=c(-FC,FC),linetype=2)+ # 
  theme(legend.background=element_rect(fill="transparent"),axis.line = element_line(color = "black"),panel.background=element_rect(fill="transparent"))#+


###Figure 2B
dataset="Melanoma-GSE91061"
infile=paste(dataset,"_diff.txt",sep="")
data5=read.table(infile,sep="\t")
pvalue=0.05
library(dplyr)
top_up5=top_n(data5[which(data5[,4]< pvalue),],n=num,wt=logFC)
top_d5=top_n(data5[which(data5[,4]< pvalue),],n=-num,wt=logFC)
library(fgsea)
library(msigdb)
library(msigdbr)
library("clusterProfiler")
m_t2g_h=c()
m_t2g_h <- msigdbr(species = "Homo sapiens", category = c("H")) %>% dplyr::select(gs_name, gene_symbol)
#m_t2g_h <- m_t2g_h%>% mutate(gs_name = paste0("H_", gs_name))
m_t2g_h=as.matrix(m_t2g_h)
library("GSVA")
em=enricher(c(rownames(top_up5),rownames(top_d5)), TERM2GENE=m_t2g_h)
temp=as.data.frame(em)
temp$function.gene.num = unlist(strsplit(temp[,"BgRatio"],"/"))[seq(1,length(unlist(strsplit(temp[,"BgRatio"],"/"))),2)]
temp$GeneRatio = temp$Count / as.numeric(temp$function.gene.num)
temp=temp[order(temp$GeneRatio),]
temp$Description=factor(temp$Description,levels=(temp$Description))
library("ggplot2")
gap = (max(temp$GeneRatio) - min(temp$GeneRatio))/5
temp$Description=gsub("HALLMARK_","",temp$Description)
p=ggplot(temp,mapping = aes(x=GeneRatio,y=Description))+
geom_point(aes(size=Count,color=p.adjust))+theme_bw()+
#scale_color_gradient(low = "blue", high = "red")+
scale_fill_viridis_c()+
scale_x_continuous(expand = c(gap, gap))+ylab(NULL)+
scale_y_discrete(limits = temp$Description[order(temp$GeneRatio)])



###########Figure 2C-D
data=read.table("Melanoma-GSE91061_diff.txt",sep="\t")
library(dplyr)
pvalue=0.05
num=200
top_up=top_n(data[which(data[,4]< pvalue),],n=num,wt=logFC)
top_d=top_n(data[which(data[,4]< pvalue),],n=-num,wt=logFC)

data1=read.table("Melanoma-GSE115821_diff.txt",sep="\t")
top_up1=top_n(data1[which(data1[,4]< pvalue),],n=num,wt=logFC)
top_d1=top_n(data1[which(data1[,4]< pvalue),],n=-num,wt=logFC)
data2=read.table("Melanoma-GSE100797_diff.txt",sep="\t")
top_up2=top_n(data2[which(data2[,4]< pvalue),],n=num,wt=logFC)
top_d2=top_n(data2[which(data2[,4]< pvalue),],n=-num,wt=logFC)

data3=read.table("Melanoma-PRJEB23709_diff.txt",sep="\t")
data4=read.table("STAD-PRJEB25780_diff.txt",sep="\t")
data5=read.table("NSCLC_GSE135222_diff.txt",sep="\t")
top_up3=top_n(data3[which(data3[,4]< pvalue),],n=num,wt=logFC)
top_d3=top_n(data3[which(data3[,4]< pvalue),],n=-num,wt=logFC)
top_up4=top_n(data4[which(data4[,4]< pvalue),],n=num,wt=logFC)
top_d4=top_n(data4[which(data4[,4]< pvalue),],n=-num,wt=logFC)
top_up5=top_n(data5[which(data5[,4]< pvalue),],n=num,wt=logFC)
top_d5=top_n(data5[which(data5[,4]< pvalue),],n=-num,wt=logFC)
gene=c(rownames(top_up),rownames(top_d))
gene1=c(rownames(top_up1),rownames(top_d1))
gene2=c(rownames(top_up2),rownames(top_d2))
gene3=c(rownames(top_up3),rownames(top_d3))
gene4=c(rownames(top_up4),rownames(top_d4))
gene5=c(rownames(top_up5),rownames(top_d5))
jarccard<-function(ge1,ge2){
	n1=length(intersect(ge1,ge2))
	n2=n1/length(unique(c(ge1,ge2)))
	return(n2)
}
alllist=list()
alllist[[1]]=gene
alllist[[2]]=gene1
alllist[[3]]=gene2
alllist[[4]]=gene3
alllist[[5]]=gene4
alllist[[6]]=gene5

out=c()
for(i in 1:6){
	for(j in 1:6){
		out=c(out,jarccard(alllist[[i]],alllist[[j]]))
	}
}
outmatrix=matrix(out,nrow=6)
#length(intersect(rownames(top_up),rownames(top_up2)))
#length(intersect(rownames(top_d),rownames(top_d2)))

library("clusterProfiler")
library(msigdbr)
library(msigdb)
library(fgsea)
m_t2g_h <- msigdbr(species = "Homo sapiens", category = c("H")) %>%dplyr::select(gs_name, gene_symbol)
m_t2g_h <- m_t2g_h%>% mutate(gs_name = paste0("H_", gs_name))
m_t2g_h=as.matrix(m_t2g_h)
em=enricher(c(rownames(top_up),rownames(top_d)), TERM2GENE=m_t2g_h)
em1=enricher(c(rownames(top_up1),rownames(top_d1)), TERM2GENE=m_t2g_h)
em2=enricher(c(rownames(top_up2),rownames(top_d2)), TERM2GENE=m_t2g_h)
em3=enricher(c(rownames(top_up3),rownames(top_d3)), TERM2GENE=m_t2g_h)
em4=enricher(c(rownames(top_up4),rownames(top_d4)), TERM2GENE=m_t2g_h)
em5=enricher(c(rownames(top_up5),rownames(top_d5)), TERM2GENE=m_t2g_h)

alllist1=list()
alllist1[[1]]=em$ID
alllist1[[2]]=em1$ID
alllist1[[3]]=em2$ID
alllist1[[4]]=em3$ID
alllist1[[5]]=em4$ID
alllist1[[6]]=em5$ID
out1=c()
for(i in 1:6){
	for(j in 1:6){
		out1=c(out1,jarccard(alllist1[[i]],alllist1[[j]]))
	}
}
out1matrix=matrix(out1,nrow=6)

library(corrplot)
diag(out1matrix) <- 0
diag(outmatrix) <- 0.4
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE","#4393C3", "#2166AC",
                           "#053061"))(100) %>% rev()
corrplot(out1matrix, is.corr = FALSE,
         method="circle", type='upper', outline = FALSE,
         tl.pos = 'd', tl.col = "black", tl.offset = 1, tl.cex = 1,
         cl.pos = 'r',
         addCoef.col = T,
         order = "original",
         number.font = 6, col = col2)
dev.off()


corrplot(outmatrix, is.corr = FALSE,
         method="circle", type='upper', outline = FALSE,
         tl.pos = 'd', tl.col = "black", tl.offset = 1, tl.cex = 1,
         cl.pos = 'r',
         addCoef.col = T,
         order = "original",
         number.font = 6, col = col2)
dev.off()





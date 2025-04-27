#####Figure 7A
pheatmap(as.matrix(plotdata),annotation_col = colgroup,annotation_colors=ann_colors,annotation_row=rowgroup,
border_color="white",cellwidth=5,cellheight=5,
cluster_rows =F,cluster_cols =T,fontsize=6,show_colnames=F,breaks=c(bk,bk1),
color=c(c2,c1)
)
########Figure 7E
bioForest=function(rt,forestFile,forestCol="red"){
        #读取输入文件
        hr <- sprintf("%.3f",rt[,1])
        hrLow  <- sprintf("%.3f",rt[,3])
        hrHigh <- sprintf("%.3f",rt[,4])
        Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
        pVal <- ifelse(rt[,2]<0.001, "<0.001", sprintf("%.3f", rt[,2]))
        gene=rownames(rt)
        #输出图形
        pdf(file=forestFile, width=6, height=4.3)
        n <- nrow(rt)
        nRow <- n+1
        ylim <- c(1,nRow)
        layout(matrix(c(1,2),nc=2),width=c(3,2.5))
        
        #绘制森林图左边的临床信息
        xlim = c(0,3)
        par(mar=c(4,2.5,2,1))
        plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
        text.cex=0.5
        text(0,n:1,gene,adj=0,cex=text.cex)
        text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);
		text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
        text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);
		text(3.1,n+1,'logOR(95% CI)',cex=text.cex,font=2,adj=1)
        
        #绘制森林图
        par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
        xlim = c(min(as.numeric(hrLow),as.numeric(hrHigh)),max(as.numeric(hrLow),as.numeric(hrHigh)))
        plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="LogOR(95% CI)")
        arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
        abline(v=0,col="black",lty=2,lwd=2)
        boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
        points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
        axis(1)
        dev.off()
}

bioForest(use,"VCAM1_forest.pdf")
#####Figure 7F
library(clusterProfiler)
library(enrichplot)
out=GSEA(fc,TERM2GENE=gmt_data,verbose = T,minGSSize=1)
gseaplot2(out,geneSetID = rownames(g)[11],
          title = "",#标题
          color = "#0070C0",#颜色
          base_size = 16,#基础大小
          rel_heights = c(0.6, 0.1, 0.3),#小图相对高度
          subplots = 1:2,#展示小图
          pvalue_table = F,#p值表格
          ES_geom = "line"#line or dot
)



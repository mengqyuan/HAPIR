####Figure 6A
library("ggplot2")
ggplot(df,aes(x=Predicted,y=T.Cells.CD8))+
  geom_smooth(method="lm",  
              se=T,         
              color="#AEC5B5")+
  geom_point(colour="#0D5122",fill="#AEC5B5",shape=21,size=3)+             
  #scale_color_manual(values=c("red"))+scale_fill_manual(values=c("#85A484"))+
  annotate(
    "text", x = min(df,na.rm = T), y = max(df,na.rm = T)-0.5,                   
    label = paste0("R=",round(-0.2848,2),                                   
                   ",P=",format(0.00354, scientific = TRUE, digits = 4)),  
    size = 5,hjust=0,vjust=1)+
  theme_classic()


####Figure 6B-D
library("ggplot2")

p=ggplot(all,mapping = aes(x=pvalue,y=gene))+
geom_point(aes(size=1,color=cor))+theme_bw()+
#scale_color_gradient(low = "blue", high = "red")+
#scale_fill_viridis_c()+
scale_colour_gradient2(low = "#002D63", mid = "#A6A6A6", high = "#A21017")+
#scale_x_continuous(expand = c(gap, gap))+
ylab(NULL)+xlab("-log10(pvalue)")+
scale_y_discrete(limits = all[,1])

score_gene_integrate<-function(gene,i,method,up)
{
	
	
	exp=matrix(nrow=56269,ncol=2)
	infile=paste(i,".Response.tsv",sep="")
	data=read.table(infile,sep="\t")
	rownames(data)=data[,1]
	data=data[,-1]
	exp=cbind(exp,data)
	exp=exp[,-c(1,2)]
	exp[is.na(exp)]=0
	exp=log2(exp+1)
	
	non_response=c()
	response=c()
	num_stat=c()

	infile=paste(i,".Response.tsv",sep="")
	data=read.table(infile,sep="\t",comment.char ="@")
	non=data[which(data[,"response_NR"]=="N"),1]
	res=data[which(data[,"response_NR"]=="R"),1]
	

	non_response=c(non_response,non)
	response=c(response,res)
	
	res_matrix1=exp[,intersect(non_response,colnames(exp))]
	sen_matrix1=exp[,intersect(response,colnames(exp))]
	train_matrix=cbind(res_matrix1,sen_matrix1)
	label=c(rep(1,dim(res_matrix1)[2]),rep(0,dim(sen_matrix1)[2]))
	if(method=="meanscore"){
		train_matrix1=t(apply(train_matrix, 1, scale))
		train_matrix1[is.na(train_matrix1)]=0
		colnames(train_matrix1)=colnames(train_matrix)
		train_matrix1=as.data.frame(train_matrix1)
		out_score=c()
		for(j in 1:length(gene)){
			genelist=unlist(strsplit(gene[j],"/"))
			up_gene=intersect(genelist,up)
			down_gene=setdiff(genelist,up)
			up_temp=train_matrix1[intersect(up_gene,rownames(train_matrix1)),]
			down_temp=train_matrix1[intersect(down_gene,rownames(train_matrix1)),]
			score=(colSums(up_temp)-colSums(down_temp))/length(genelist)
			score[is.na(score)]=0
			out_score=cbind(out_score,score)
		}
		colnames(out_score)=names(gene)
		rownames(out_score)=colnames(train_matrix1)
        }
	if(method=="gsva"){
		library("GSVA")
		geneset=list()
		for(j in 1:length(gene)){
			genelist=as.vector(unlist(strsplit(gene[j],"/")))
			geneset[[names(gene)[j]]]=genelist
		}
		re=gsva(as.matrix(train_matrix),geneset,method = "gsva",min.sz= 1,max.sz=1000,ssgsea.norm= T,parallel.sz= 1L)
		temp=setdiff(names(geneset),rownames(re))
		buquan=matrix(0,nrow=length(temp),ncol=dim(re)[2])
		rownames(buquan)=temp
		out_score=t(rbind(re,buquan))
	}
	if(method=="ssgsea"){
		library("GSVA")
		geneset=list()
		for(j in 1:length(gene)){
			genelist=as.vector(unlist(strsplit(gene[j],"/")))
			geneset[[names(gene)[j]]]=genelist
		}
		re=gsva(as.matrix(train_matrix),geneset,method = "ssgsea",min.sz= 1,max.sz=1000,ssgsea.norm= T,parallel.sz= 1L)
		temp=setdiff(names(geneset),rownames(re))
		buquan=matrix(0,nrow=length(temp),ncol=dim(re)[2])
		rownames(buquan)=temp
		out_score=t(rbind(re,buquan))
	}
	if(method=="aucell"){
		library("AUCell")
		geneset=list()
		for(j in 1:length(gene)){
			genelist=as.vector(unlist(strsplit(gene[j],"/")))
			geneset[[names(gene)[j]]]=genelist
		}
		cells_rankings <- AUCell_buildRankings(as.matrix(train_matrix))
		cells_AUC <- AUCell_calcAUC(geneset, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
		aucs <- getAUC(cells_AUC)
		temp=setdiff(names(geneset),rownames(aucs))
		buquan=matrix(0,nrow=length(temp),ncol=dim(aucs)[2])
		rownames(buquan)=temp
		out_score=t(rbind(aucs,buquan))
	}
	if(method=="addmodulescore"){
		library("Seurat")
		ob=CreateSeuratObject(counts=as.matrix(train_matrix),min.cells = 0,min.features = 0)
		for(j in 1:length(gene)){
			genelist=list(as.vector(unlist(strsplit(gene[j],"/"))))
			if(length(intersect(rownames(ob),genelist[[1]]))>0){
			ob=AddModuleScore(ob,features=genelist,slot ="count",name = names(gene)[j])
			}
		}
		out_score=ob@meta.data[,4:dim(ob@meta.data)[2]]
		newname=paste(names(gene),"1",sep="")
		temp=setdiff(newname,colnames(out_score))
		if(length(temp)>0){
		buquan=matrix(0,ncol=length(temp),nrow=dim(out_score)[1])
		colnames(buquan)=temp
		out_score=cbind(out_score,buquan)}
	}

	return(list(score=out_score,label=label))
}


external_validation=function(model,test,label,method)
{
		if(method=="logit"){
		pred <- predict(model, newdata = as.data.frame(test), type = "response")
	}
	
	if(method=="Enet"){
		pred <- predict(model,as.data.frame(test),type="prob")
		pred=pred[,2]
	}
	
	if(method=="Ridge"){
		pred <- predict(model,test,type="response")
	}
	
	if(method=="RF"){
		pred <- predict(model,as.data.frame(test),type="prob")
		pred=pred[,2]
	}
	
	if(method=="xgboost"){
		dtest <- xgb.DMatrix(data=test,label=label)
		pred <- predict(model, newdata=dtest,type="response")
	}
	
	if(method=="SVM"){
		pred1 <- predict(model, as.data.frame(test),probability =TRUE)
		pred <- attr(pred1,"probabilities")[,1]
	}
	if(method=="knn"){
		pred1 <- predict(model,as.data.frame(test),type="prob")
		pred <- pred1[,2]
	}
        rocobj <- roc(as.numeric(label),pred,direction="<")
        auc <- round(auc(rocobj),4)
        return(auc)
}


single_train_model<-function(dataset,num=200,method,score_method,sig){
infile=paste(dataset,"_diff.txt",sep="")
data5=read.table(infile,sep="\t") ###
pvalue=0.05
library(dplyr)
top_up5=top_n(data5[which(data5[,4]< pvalue),],n=num,wt=logFC)
top_d5=top_n(data5[which(data5[,4]< pvalue),],n=-num,wt=logFC)
library(fgsea)
library(msigdb)
library(msigdbr)
library("clusterProfiler")
m_t2g_h=c()
if(sig=="hallmark"){
m_t2g_h <- msigdbr(species = "Homo sapiens", category = c("H")) %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g_h <- m_t2g_h%>% mutate(gs_name = paste0("H_", gs_name))
m_t2g_h=as.matrix(m_t2g_h)
}

library("GSVA")
em=enricher(c(rownames(top_up5),rownames(top_d5)), TERM2GENE=m_t2g_h)
all_list=em$geneID
temp=em$ID
temp=gsub(".","",temp,fixed=T)
temp=gsub(")","",temp,fixed=T)
temp=gsub("(","",temp,fixed=T)
temp=gsub(" ","_",temp)
temp=gsub("-","_",temp)
names(all_list)=temp

score=score_gene_integrate(all_list,dataset,score_method,rownames(top_up5))
library(glmnet)
library(caret)
library(randomForest)
library(xgboost)
library(e1071)

traindata=cbind(score$label,score$score)
colnames(traindata)[1]="label"

if(method=="logit"){
	model<- glm(label~.,family=binomial(link = "logit"),data = as.data.frame(traindata))
}


if(method=="Enet"){
	traindata=as.data.frame(traindata)
	traindata$label <- factor(traindata$label)
	model <- train(
		 label~.,data = traindata, method="glmnet", 
		 trControl = trainControl("cv",number=10),
		 tuneLength = 10)
}

if(method=="Ridge"){
	glmnet1 <- cv.glmnet(x=as.matrix(traindata[,-1]),y=as.character(traindata[,1]),alpha=0,family="binomial")
	model <- glmnet(x=traindata[,-1],y=as.character(traindata[,1]),alpha=0,lambda = glmnet1$lambda.min,family="binomial")
}

if(method=="RF"){
	traindata=as.data.frame(traindata)
	traindata$label <- factor(traindata$label)
	model <- randomForest(label~., data = traindata, importance = TRUE)
}

if(method=="xgboost"){
	dtrain <- xgb.DMatrix(data=as.matrix(traindata[,-1]),label=as.character(traindata[,1]))
	model <- xgboost(data=dtrain,nround=100,objective="binary:logistic",learning_rate=0.1)
}

if(method=="SVM"){
	traindata=as.data.frame(traindata)
	traindata$label <- factor(traindata$label)
	model <- svm(label ~. , traindata,probability = TRUE)
}
if(method=="knn"){
	train_1 <- apply(traindata[,-1],2,scale)
	traindata=as.data.frame(cbind(traindata[,1],train_1))
	colnames(traindata)[1]="label"
	traindata$label <- factor(traindata$label)
	grid1 <- expand.grid(.k = seq(2, 20, by = 1))
	model <- train(
		label~.,data = traindata, method="knn", 
		trControl = trainControl("cv",number=10),
		tuneGrid = grid1)
}
library(pROC)
library(ggplot2)

all_dataset=c("NSCLC_GSE135222")

all_auc=c()
for(i in sort(all_dataset))
{
	print(i)
	datasetexp=score_gene_integrate(all_list,i,score_method,rownames(top_up5))
	if(method=="xgboost"){
	auc=external_validation(model,as.matrix(datasetexp$score[,model$feature_names]),datasetexp$label,method) 
	}
	if(method!="xgboost"){
	auc=external_validation(model,as.matrix(datasetexp$score),datasetexp$label,method) 
	}
	all_auc=rbind(all_auc,c(i,as.matrix(table(datasetexp$label))[2],as.matrix(table(datasetexp$label))[1],auc))
}

colnames(all_auc)=c("dataset","NR","R","auc")
}


genenum=c(100,200,300,400,500)
model=c("logit","Enet","Ridge","RF","xgboost","SVM","knn")
score=c("gsva","ssgsea","aucell","addmodulescore","meanscore")
for(i in genenum){
	for(j in model){
		for(m in score){
			for(k in geneset){
				single_train_model("Melanoma-GSE91061",i,j,m,k)
			}
		}
	}
}

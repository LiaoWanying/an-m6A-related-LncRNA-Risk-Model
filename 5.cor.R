#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
corFilter=0.4            #相关系数过滤标准
pvalueFilter=0.001       #p值过滤标准
setwd("C:\\m6AlncRNA\\09.m6aLncExp")     #设置工作目录

#读取lncRNA表达文件,并对数据进行处理
rt=read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
sampleType=c(rep(1,conNum), rep(2,treatNum))

#读取m6a基因表达文件,并对数据进行处理
rt1=read.table("m6aGeneExp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
m6A=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
m6A=avereps(m6A)
m6A=m6A[rowMeans(m6A)>0.1,]

#删掉正常样品
group=sapply(strsplit(colnames(m6A),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
m6A=m6A[,group==0]

#相关性检验
outTab=data.frame()
for(i in row.names(lncRNA)){
	if(sd(lncRNA[i,])>0.1){
		test=wilcox.test(data[i,] ~ sampleType)
		if(test$p.value<0.05){
			for(j in row.names(m6A)){
				x=as.numeric(lncRNA[i,])
				y=as.numeric(m6A[j,])
				corT=cor.test(x,y)
				cor=corT$estimate
				pvalue=corT$p.value
				if((cor>corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(m6A=j,lncRNA=i,cor,pvalue,Regulation="postive"))
				}
				if((cor< -corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(m6A=j,lncRNA=i,cor,pvalue,Regulation="negative"))
				}
			}
		}
	}
}

#输出相关性网络
write.table(file="net.network.txt",outTab,sep="\t",quote=F,row.names=F)
#输出相关性节点属性
lncNode=data.frame(Node=unique(as.vector(outTab[,"lncRNA"])), Type="lncRNA")
mrnaNode=data.frame(Node=unique(as.vector(outTab[,"m6A"])), Type="m6A")
nodeOut=rbind(lncNode, mrnaNode)
write.table(nodeOut, file="net.node.txt", sep="\t", quote=F, row.names=F)

#输出m6a相关lncRNA表达量
m6aLncRNA=unique(as.vector(outTab[,"lncRNA"]))
m6aLncRNAexp=data[m6aLncRNA,]
m6aLncRNAexp=rbind(ID=colnames(m6aLncRNAexp), m6aLncRNAexp)
write.table(m6aLncRNAexp,file="m6aLncExp.txt",sep="\t",quote=F,col.names=F)




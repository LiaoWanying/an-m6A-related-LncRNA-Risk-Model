#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)            #引用包
expFile="symbol.txt"      #表达输入文件
geneFile="gene.txt"       #基因列表文件
setwd("C:\\m6AlncRNA\\08.m6aExp")    #设置工作目录

#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#获取m6a基因的表达量
gene=read.table(geneFile, header=T, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#输出结果
out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="m6aGeneExp.txt",sep="\t",quote=F,col.names=F)
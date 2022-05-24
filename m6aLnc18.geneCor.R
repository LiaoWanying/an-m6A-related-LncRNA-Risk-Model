
library(limma)
library(corrplot)
setwd("/Users/lwy/Desktop/重跑v2/geneCor")
expFile="symbol.txt"       
lncFile="uniSigExp.txt"     
geneFile="gene.txt"


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]


gene=read.table(geneFile, header=T, sep="\t", check.names=F, row.names=1)
samemRNA=intersect(rownames(gene),rownames(data))
datamRNA=data[samemRNA,]


lncRT=read.table(lncFile, header=T, sep="\t", check.names=F, row.names=1)
lncRNA=colnames(lncRT)[3:ncol(lncRT)]
sameGene=intersect(lncRNA, row.names(data))
datalncRNA=data[sameGene,]


datamRNA=t(datamRNA)
datalncRNA=t(datalncRNA)
M=cor(datamRNA,datalncRNA,method="pearson")
p.mat <- cor.mtest(data)
p1=p.mat$p[1:23,24:52]


pdf(file="cor.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "full",
         tl.cex=0.8, pch=T,
         p.mat = p1,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")


#???ð?
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
lncFile="uniSigExp.txt"     #Ԥ?????ص?lncRNA?б?
expFile="lncRNA.txt"        #lncRNA?????ļ?
setwd("/Users/lwy/Desktop/重跑/13")       #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#??ȡԤ??lncRNA????��
lncRNA=read.table(lncFile, header=T, sep="\t", check.names=F, row.names=1)
data=data[colnames(lncRNA)[3:ncol(lncRNA)],]
exp=data

#????????????Ŀ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       #????????Ʒ??Ŀ
treatNum=length(group[group==0])     #????????Ʒ??Ŀ
sampleType=c(rep(1,conNum), rep(2,treatNum))

#????????????
sigVec=c()
for(i in row.names(data)){
	test=wilcox.test(data[i,] ~ sampleType)
	pvalue=test$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(i, Sig))
}
row.names(data)=sigVec

#??ͼ???ӻ?
Type=c(rep("Normal",conNum), rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
data=log2(data+1)
pdf("heatmap.pdf", width=7.5, height=4.7)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

#??????ת????ggplot2?????ļ?
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#????????ͼ
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
	     ylab="Gene expression",
	     xlab="",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#????????ͼ
pdf(file="boxplot.pdf", width=7.5, height=5)
print(p1)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

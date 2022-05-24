#install.packages("ggpubr")


library(ggpubr)                    #???ð?
tciaFile="TCIA.txt"                #???????ƴ????ļ?
scoreFile="allRisk.txt"     #m6A???ַ????ļ?
setwd("/Users/lwy/Desktop/重跑v2/data16")     #?޸Ĺ???Ŀ¼

#??ȡ???????ƴ????ļ?
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#??ȡm6A???ַ????ļ?
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
sameSample=intersect(row.names(ips), row.names(score))
ips=ips[sameSample, , drop=F]
score=score[sameSample, "group", drop=F]
data=cbind(ips, score)

#???ñȽ???
data$group=factor(data$group, levels=c("low", "high"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#?????????ƴ??ֽ???ѭ??,?ֱ?????С????ͼ
for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "group")]
	colnames(rt)=c("IPS", "group")
	gg1=ggviolin(rt, x="group", y="IPS", fill = "group", 
	         xlab="m6Ascore", ylab=i,
	         legend.title="m6Ascore",
	         palette=c("#0066FF", "#FF0000"),
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0(i, ".pdf"), width=6, height=5)
	print(gg1)
	dev.off()
}


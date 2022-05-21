######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("pheatmap")


library(pheatmap)         #???ð?
setwd("/Users/lwy/Desktop/重跑/29")      #???ù???Ŀ¼

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
	rt=rt[order(rt$riskScore),]      #????riskScore????Ʒ????
		
	#???Ʒ???????
	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	line[line>10]=10
	pdf(file=riskScoreFile, width=7, height=4)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)", ylab="Risk score",
		 col=c(rep("green",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
	dev.off()
		
	#????????״̬ͼ
	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="green"
	pdf(file=survStatFile, width=7, height=4)
	plot(rt$futime, pch=19,
		 xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
		
	#???Ʒ?????ͼ
	rt1=rt[c(3:(ncol(rt)-2))]
	rt1=t(rt1)
	annotation=data.frame(type=rt[,ncol(rt)])
	rownames(annotation)=rownames(rt)
	pdf(file=heatmapFile, width=7, height=4)
	pheatmap(rt1, 
		     annotation=annotation, 
		     cluster_cols = FALSE,
		     cluster_rows = FALSE,
		     show_colnames = F,
		     scale="row",
		     color = colorRampPalette(c(rep("green",3.5), "white", rep("red",3.5)))(50),
		     fontsize_col=3,
		     fontsize=7,
		     fontsize_row=8)
	dev.off()
}
#tarin??????????
bioRiskPlot(inputFile="allRisk.txt",riskScoreFile="all.riskScore.pdf",survStatFile="all.survStat.pdf",heatmapFile="all.heatmap.pdf")
#test??????????
bioRiskPlot(inputFile="testRisk.txt",riskScoreFile="test.riskScore.pdf",survStatFile="test.survStat.pdf",heatmapFile="test.heatmap.pdf")
bioRiskPlot(inputFile="trainRisk.txt",riskScoreFile="train.riskScore.pdf",survStatFile="train.survStat.pdf",heatmapFile="train.heatmap.pdf")


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

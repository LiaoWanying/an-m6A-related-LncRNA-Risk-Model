
#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


#???ð?
library(survival)
library(survminer)
library(timeROC)
riskFile="testRisk.txt"         #?????????ļ?
cliFile="clinical.txt"      #?ٴ??????ļ?
setwd("/Users/lwy/Desktop/重跑/ROC")     #?޸Ĺ???Ŀ¼

#??ȡ?????????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#??ȡ?ٴ??????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#??????ɫ
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######????1 2 3????ROC????######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
	           marker=risk$riskScore,cause=1,
	           weighting='aalen',
	           times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
	   c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	     paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	     paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	   col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()


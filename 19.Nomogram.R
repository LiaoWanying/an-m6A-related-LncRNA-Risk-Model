install.packages("rms")

library(rms)
setwd("/Users/lwy/Desktop/重跑v2/data13/1")                         #???ù???Ŀ¼

#TCGA????ͼ????
riskFile="allRisk.txt"
cliFile="clinical.txt"
outFile="Nomogram.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        #??ȡ?????ļ?
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          #??ȡ?ٴ??ļ?
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
#???ݴ???
dd <- datadist(rt)
options(datadist="dd")
#???ɺ???
f <- cph(Surv(futime, fustat) ~stage+riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
#建立nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x), function(x) surv(4, x), function(x) surv(5, x)), 
                lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival", "4-year survival", "5-year survival"), 
                maxscale=100, 
                fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
#nomogram可视化
pdf(file=outFile,height=7.5,width=14)
plot(nom)
dev.off()


install.packages("survival")
install.packages("survminer")


#???ð?
library(survival)
library(survminer)
riskFile="allRisk.txt"      #?????????ļ?
cliFile="clinical.txt"      #?ٴ??????ļ?
setwd("/Users/lwy/Desktop/重跑/data10")                   #???ù???Ŀ¼
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #??ȡ?ٴ??ļ?

#???ݺϲ?
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
data=cbind(futime=risk[,1],fustat=risk[,2],cli,risk=risk[,"risk"])

#???ٴ???Ϣ????ѭ??
for(i in colnames(data[,3:(ncol(data)-1)])){
    rt=data[,c("futime","fustat",i,"risk")]
    rt=rt[(rt[,i]!="unknow"),]
    colnames(rt)=c("futime","fustat","clinical","risk")
	tab=table(rt[,"clinical"])
	tab=tab[tab!=0]
	#??ÿ???ٴ???Ϣ?????ķ???????ѭ??
	for(j in names(tab)){
		rt1=rt[(rt[,"clinical"]==j),]
		tab1=table(rt1[,"risk"])
		tab1=tab1[tab1!=0]
		labels=names(tab1)
		if(length(labels)==2){
			titleName=j
			if((i=="age") | (i=="Age") | (i=="AGE")){
				titleName=paste0("age",j)
			}
			diff=survdiff(Surv(futime, fustat) ~risk,data = rt1)
			pValue=1-pchisq(diff$chisq,df=1)
			if(pValue<0.001){
				pValue="p<0.001"
			}else{
				pValue=paste0("p=",sprintf("%.03f",pValue))
			}
			fit <- survfit(Surv(futime, fustat) ~ risk, data = rt1)
			#????????????
			surPlot=ggsurvplot(fit, 
			           data=rt1,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           title=paste0("Patients with ",titleName),
			           legend.title="Risk",
			           legend.labs=labels,
			           font.legend=12,
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette=c("red", "blue"),
			           risk.table=TRUE,
			       	   risk.table.title="",
			           risk.table.col = "strata",
			           risk.table.height=.25)
		    #????ͼƬ
		    j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
			pdf(file=paste0("survival.",i,"_",j,".pdf"),onefile = FALSE,
			       width = 6,        #ͼƬ?Ŀ???
			       height =5)        #ͼƬ?ĸ߶?
			print(surPlot)
			dev.off()
		}
	}
}


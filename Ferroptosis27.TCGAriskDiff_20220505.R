
library(limma)

expFile="mRNA.txt"        
riskFile="allRisk.txt"      
fdrFilter=0.05         
logFCfilter=1             
setwd("D:\\biowolf\\Ferroptosis\\27.riskDiff\\TCGA")  
#setwd("C:/data/expl/order/20220505793/")


rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data=avereps(t(data))
data=t(data)


risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
sameSample=intersect(colnames(data),row.names(risk))
data=data[,sameSample]
risk=risk[sameSample,]


riskLow=risk[risk$risk=="low",]
riskHigh=risk[risk$risk=="high",]
dataLow=data[,row.names(riskLow)]
dataHigh=data[,row.names(riskHigh)]
data=cbind(dataLow,dataHigh)
data=data[rowMeans(data)>1,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep(1,conNum),rep(2,treatNum))

# 导入临床数据
risk[1:6,]
data[1:6, 1:6]
inFile <- "clinicalpractice.xls"
samp=read.table(inFile,header=T,sep="\t",row.names=1,check.names=F)
samp[1:6,]

# high和low差异分析
# Limma差异表达分析
library(limma)
eset <- data
eset[1:6, 1:6]
clindat <- cbind(risk[colnames(eset),], samp[colnames(eset),])
clindat[1:6,]
table(clindat[, "risk"])
design <- model.matrix(~ 0 + risk + gender + age, data=clindat)
colnames(design)
colnames(design) <- c("high", "low", "male", "age")
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(
    high-low,
    levels=design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2
results <- decideTests(fit2)
vennDiagram(results)
outTab <- topTable(fit2, coef=1, adjust="BH", number=nrow(eset))
nrow(outTab)
write.csv(outTab, "result_limma_high-low.csv")

outTab[1:6,]
outDiff=outTab[abs(outTab$logFC)>logFCfilter & outTab$adj.P.Val<fdrFilter,]
nrow(outDiff)
write.table(outDiff,file="riskDiff.txt",sep="\t",row.names=T,quote=F)

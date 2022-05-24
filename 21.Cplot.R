#Start
#options(stringsAsFactors=F)#有意思，全局设置不转换为因子
#清除所有变量
rm(list=ls())
#运行开始
timestart<-Sys.time() #记录开始运行时间
options(stringsAsFactors=F, quote="",unzip = "internal")#有意思，全局设置不转换为因子
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
workDir="/Users/lwy/Desktop/重跑v2/data13/3"        #工作目录
setwd(workDir)
list.files() #查看工作目录下的文件
File1_0<-read.table("expTime.txt",header=T,row.names=1, sep="\t",check.names = F, quote="")#quote=""
str(File1_0)
# 'data.frame':   163 obs. of  4 variables:
 # $ futime   : num  0.00274 0.03288 0.03288 0.03836 0.04384 ...
 # $ fustat   : int  0 0 1 1 1 0 1 0 1 1 ...
 # $ RiskScore: num  0.107 0.538 1.934 -3.284 -0.141 ...
 # $ Stage    : int  0 0 0 0 0 1 0 0 0 0 ...
#得转化为因子
File1_0$Stage <- factor (File1_0$Stage)
str(File1_0)
abc <- File1_0



# ###goodness of fit_______________####这一段不知道有何用处，出了一堆奇奇怪怪的散点图
# library(survival)
# library(survminer)
# fit<-coxph(Surv(futime,fustat)~RiskScore+Stage , data = abc)
# ggcoxdiagnostics(fit, type = "schoenfeld")
# ggcoxdiagnostics(fit, type = "schoenfeld", ox.scale="time")
# ggcoxdiagnostics(fit, type = "deviance")
# ggcoxdiagnostics(fit, type = "martingale")

####C-index#####
abc <- File1_0
names(abc)
library(rms)
library(survival)
set.seed(1)# the method of validate is using random
fit.cph <- cph(Surv(futime, fustat)~ RiskScore+Stage,data=abc, 
               x=TRUE,y=TRUE,surv=TRUE)

fit.cph

v<-validate(fit.cph, dxy=TRUE, B=1000)

Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"]

bias_corrected_c_index  <- abs(Dxy)/2+0.5
orig_c_index <- abs(orig_Dxy)/2+0.5
bias_corrected_c_index    #HCC数据显示为 0.6863544
orig_c_index              #HCC数据显示为 0.6935687

library(rms)
#install.packages("CsChange")
library(CsChange)    ###how c-index changed from model1 to model2###
dd<- datadist(abc)
options(datadist="dd")
model1=cph(Surv(futime, fustat)~ RiskScore+Stage,abc)
model2=cph(Surv(futime, fustat)~ Stage,abc)
CsChange(model1,model2,data=abc,nb=500)
#这一段也很奇怪

source("stdca.R",encoding="utf-8") 
####nomogram__________________________######
library(nomogramEx)
library(rms)
abc <- File1_0
attach(abc)
ddist<-datadist(Stage, RiskScore )
options(datadist='ddist')
f<- cph(formula(Surv(futime, fustat)~Stage+RiskScore),#此处决定列线图各变量顺序
		data = abc, x=TRUE, y=TRUE, surv=TRUE, time.inc=3) #哪怕前面有了attach，还是需要有data=abc数据框
		#data = abc, 
surv <- Survival(f)

nomo <- nomogram(f, fun=list(function(x) surv(1,x),function(x) surv(3,x),function(x) surv(5,x)),
                 lp=FALSE,funlabel=c("1-year Survival Prob","3-year Survival Prob","5-year Survival Prob"))
str(nomo)
pdf(file = "LIHC-Nomo-lwd=20.pdf"
	, width = 21
	, height = 8
	#, onefile
	#, family = "Times"#默认"Helvetica"无衬线体sans-serif font，
	                #可以改为"Times"衬线体，monospaced font (to "Courier").
	# , title
	# , fonts
	# , version
	# , paper
	# , encoding
	# , bg
	# , fg
	# , pointsize
	# , pagecentre
	# , colormodel
	# , useDingbats
	# , useKerning
	# , fillOddEven
	# , compress
	)
plot(nomo,cex.axis = 1.2,cex.var = 1.2, lwd=20,
     label.every=2, 
     
     col.grid=gray(c(0.8,0.95))
     )

dev.off()
#
detach(abc)



##c-plot___________________________________________###
library(nomogramEx)
library(rms)
library(survival)


ddist <- datadist(abc)
options(datadist='ddist')

units(abc$futime) <- "Years"


#####1years________________________________####
fcox1 <- cph(Surv(futime, fustat) ~ RiskScore+Stage,
             surv=T,x=T, y=T,time.inc = 1,data=abc) #time.inc = 1表示第1年
cal1 <-  calibrate(fcox1,  cmethod="KM", method="boot", u=1, m=100, B=500)

pdf("COAD-Cplot-1y.pdf")
plot(cal1,lwd=3,lty=1,cex=0.000001,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),riskdist=F,
     xlab="Predicted Probability of 1-Year Survival",
     ylab="Actual Proportion of 1-Year Survival",cex.lab=1.5,cex.axis=1.5,
     col=c(rgb(192,98,83,maxColorValue=255)))
lines(cal1[,c("mean.predicted","KM")],type="b",lwd=3,col=c(rgb(192,98,83,maxColorValue=255)),pch=16)
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
dev.off()


#####3years________________________________####
fcox3 <- cph(Surv(futime, fustat) ~ RiskScore+Stage,
             surv=T,x=T, y=T,time.inc = 3,data=abc) 
cal3 <-  calibrate(fcox3,  cmethod="KM", method="boot", u=3, m=100, B=500)

pdf("COAD-Cplot-3y.pdf")
plot(cal3,lwd=3,lty=1,cex=0.000001,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),riskdist=F,
     xlab="Predicted Probability of 3-Year Survival",
     ylab="Actual Proportion of 3-Year Survival",cex.lab=1.5,cex.axis=1.5,
     col=c(rgb(192,98,83,maxColorValue=255)))
lines(cal3[,c("mean.predicted","KM")],type="b",lwd=3,col=c(rgb(192,98,83,maxColorValue=255)),pch=16)
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
dev.off()



#####2years________________________________####
fcox2 <- cph(Surv(futime, fustat) ~ RiskScore+Stage,
             surv=T,x=T, y=T,time.inc = 2,data=abc) 
cal2 <-  calibrate(fcox5,  cmethod="KM", method="boot", u=2, m=100, B=500)

pdf("COAD-Cplot-2y.pdf")
plot(cal5,lwd=3,lty=1,cex=0.000001,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),riskdist=F,
     xlab="Predicted Probability of 5-Year Survival",
     ylab="Actual Proportion of 5-Year Survival",cex.lab=1.5,cex.axis=1.5,
     col=c(rgb(192,98,83,maxColorValue=255)))
lines(cal5[,c("mean.predicted","KM")],type="b",lwd=3,col=c(rgb(192,98,83,maxColorValue=255)),pch=16)
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
dev.off()


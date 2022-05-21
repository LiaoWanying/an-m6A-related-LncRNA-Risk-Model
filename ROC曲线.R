setwd("/Users/lwy/Desktop/重跑v2/data13/2")
options(stringsAsFactors=FALSE)

# 2???????ݼ?????????????
bc<-read.csv("expTime.csv")
dim(bc)
#View(bc)
bc[1:5,1:5]
#??2???鿴???ݿ?C????????????
# ??????
length(bc[1, ])
# ??????
length(bc[ ,1])
# ????+??????
dim(bc)

summary(bc[,"futime"])
# ??????????ʱ??????????Ϊ??λ?ģ????Գ???365????????
bc[,"futime"]=bc[,"futime"]*365
dim(bc)
bc[1:5,1:5]
summary(bc[,"futime"])

bc <- na.omit(bc)
dim(bc)
bc[1:5:5]
#??2???鿴???ݿ?C????????????
# ??????
length(bc[1, ])
# ??????
length(bc[ ,1])
# ????+??????
dim(bc)


library(survival)
##################################################??1??ROC#######################################
# ??1????ROC
library(timeROC)
fit <- coxph(Surv(futime, fustat) ~ PathologicStage, data=bc)
# ????????
bc$pred <- predict(fit, newdata=bc)
roc3 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting='aalen',    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*2,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc3 <- roc3$AUC[2]
auc3


fit <- coxph(Surv(futime, fustat) ~ Risk, data=bc)
# Ԥ??ģ??
bc$pred <- predict(fit, newdata=bc)
roc4 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting="aalen",    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*2,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc4 <- roc4$AUC[2]
auc4

fit <- coxph(Surv(futime, fustat) ~ Risk + PathologicStage, data=bc)
# Ԥ??ģ??+????????
bc$pred <- predict(fit, newdata=bc)
roc5 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting="aalen",    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*2,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc5 <- roc5$AUC[2]
auc5

#pdf(file="ROC-1??.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
plot(roc3$FP[,2], roc3$TP[,2], type="l", xlim=c(0,1), ylim=c(0,1),col=2, 
    xlab="False positive rate", ylab="True positive rate",
    lwd=6)
lines(roc4$FP[,2], roc4$TP[,2], type="l",col=3,lwd=6)
lines(roc5$FP[,2], roc5$TP[,2], type="l",col=4,lwd=6)
abline(0,1)
legend("bottomright", legend=c(
"Pathologic Stage model",
"Prognostic model",
"Combined model"),
lty=c(1, 1, 1, 1, 1), col=c(2, 3, 4), box.col="transparent",lwd=6)
# 1??ʵ?֣?2??????
text(x=0.6, y =0.2, labels ="2 year AUC")
#dev.off()
auc3 # ????????
auc4 # Ԥ??ģ??
auc5 # Ԥ??ģ??+????????


# ??Ϊƽ??????
# https://github.com/cran/pROC/blob/master/R/smooth.R
# smooth.method="binormal", smooth.n=512
smooth.roc.binormal <- function(roc, n) {
  df <- data.frame(sp=qnorm(roc$sp * ifelse(roc$percent, 1/100, 1)), se=qnorm(roc$se * ifelse(roc$percent, 1/100, 1)))
  df <- df[apply(df, 1, function(x) all(is.finite(x))),]
  if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
  model <- lm(sp~se, df)
  if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
  se <- qnorm(seq(0, 1, 1/(n-1)))
  sp <- predict(model, data.frame(se))

  return(list(sensitivities = pnorm(se) * ifelse(roc$percent, 100, 1),
              specificities = pnorm(sp) * ifelse(roc$percent, 100, 1),
              model = model))
}

# ???????ڣ???????ģ?ͣ???Ҫ???й????????ĵ?
se1 <- roc3[["FP"]][, 2]
sp1 <- roc3[["TP"]][, 2]
se1
sp1
# ???й????????ĵ?
se11 <- seq(se1[1], se1[2], length=10)
sp11 <- (sp1[1] - sp1[2])/(se1[1] - se1[2]) * (se11 - se1[2]) + sp1[2]
se12 <- seq(se1[2], se1[3], length=10)
sp12 <- (sp1[2] - sp1[3])/(se1[2] - se1[3]) * (se12 - se1[3]) + sp1[3]
se1 <- c(se11, se12)
sp1 <- c(sp11, sp12)
#
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t1 <- se2
sp_t1 <- sp2

# Ԥ??ģ??
se1 <- roc4[["FP"]][, 2]
sp1 <- roc4[["TP"]][, 2]
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t2 <- se2
sp_t2 <- sp2

# Ԥ??ģ??+????????
se1 <- roc5[["FP"]][, 2]
sp1 <- roc5[["TP"]][, 2]
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t3 <- se2
sp_t3 <- sp2

text1 <- paste0("Pathologic Stage model: AUC = ", sprintf("%.3f", roc3[["AUC"]][2]))
text2 <- paste0("Prognostic model: AUC = ", sprintf("%.3f", roc4[["AUC"]][2]))
text3 <- paste0("Combined model: AUC = ", sprintf("%.3f", roc5[["AUC"]][2]))

outfile <- "fig_roc_smooth_year2.pdf"
pdf(outfile, onefile=FALSE, width=15/2.54, height=15/2.54)
par(mar=c(4.5, 4.5, 2, 2))
cols <- c("red", "blue", "darkgreen")
ltys <- c(1, 1, 1, 1, 1)
lwds <- c(6, 6, 6, 6, 6)
legends <- c(text1, text2, text3)
plot(
    se_t1, sp_t1, type="l", col=cols[1], lty=ltys[1], lwd=lwds[1],
    xlim=c(0, 1), ylim=c(0, 1),
    xlab="1 - Specificity",
    ylab="Sensitivity"
)
polygon(c(se_t1, 1), c(sp_t1, 0), col="lightblue")
lines(se_t2, sp_t2, col=cols[2], lty=ltys[2], lwd=lwds[2])
lines(se_t3, sp_t3, col=cols[3], lty=ltys[3], lwd=lwds[3])
abline(0, 1, col="grey")
legend(
    0.9, 0.4, legend=legends, text.col=cols, adj=1, box.col=NA
)
text(x=0.6, y=0.42, labels="1 year AUC")
dev.off()

outfile <- "fig_roc_smooth_year2.tiff"
tiff(outfile, width=15, height=15, unit="cm", res=600, compression="lzw+p")
par(mar=c(4.5, 4.5, 2, 2))
cols <- c("red", "blue", "darkgreen")
ltys <- c(1, 1, 1, 1, 1)
lwds <- c(6, 6, 6, 6, 6)
legends <- c(text1, text2, text3)
plot(
    se_t1, sp_t1, type="l", col=cols[1], lty=ltys[1], lwd=lwds[1],
    xlim=c(0, 1), ylim=c(0, 1),
    xlab="1 - Specificity",
    ylab="Sensitivity"
)
polygon(c(se_t1, 1), c(sp_t1, 0), col="lightblue")
lines(se_t2, sp_t2, col=cols[2], lty=ltys[2], lwd=lwds[2])
lines(se_t3, sp_t3, col=cols[3], lty=ltys[3], lwd=lwds[3])
abline(0, 1, col="grey")
legend(
    0.9, 0.4, legend=legends, text.col=cols, adj=1, box.col=NA
)
text(x=0.6, y=0.42, labels="2 year AUC")
dev.off()




##################################################??3??ROC#######################################
# ??3????ROC
library(timeROC)
fit <- coxph(Surv(futime, fustat) ~ PathologicStage, data=bc)
# ????????
bc$pred <- predict(fit, newdata=bc)
roc3 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting='aalen',    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*3,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc3 <- roc3$AUC[2]
auc3

fit <- coxph(Surv(futime, fustat) ~ Risk, data=bc)
# Ԥ??ģ??
bc$pred <- predict(fit, newdata=bc)
roc4 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting='aalen',    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*3,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc4 <- roc4$AUC[2]
auc4

fit <- coxph(Surv(futime, fustat) ~ Risk +PathologicStage, data=bc)
# Ԥ??ģ??+????????
bc$pred <- predict(fit, newdata=bc)
roc5 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting='aalen',    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*3,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc5 <- roc5$AUC[2]
auc5

#pdf(file="ROC-3??.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
plot(roc3$FP[, 2], roc3$TP[, 2], type="l", xlim=c(0,1), ylim=c(0,1),col=2, 
    xlab="False positive rate", ylab="True positive rate",
    lwd=6)
lines(roc4$FP[, 2], roc4$TP[, 2], type="l",col=3,lwd=6)
lines(roc5$FP[, 2], roc5$TP[, 2], type="l",col=4,lwd=6)
abline(0,1)
legend("bottomright", legend=c(
"Pathologic Stage model",
"Prognostic model",
"Combined model"),
lty=c(1, 1, 1, 1, 1), col=c(2, 3, 4), box.col="transparent",lwd=6)
# 1??ʵ?֣?2??????
text(x=0.6, y =0.2, labels ="3 year AUC")
#dev.off()
auc3 # ????????
auc4 # Ԥ??ģ??
auc5 # Ԥ??ģ??+????????


# ??Ϊƽ??????
# https://github.com/cran/pROC/blob/master/R/smooth.R
# smooth.method="binormal", smooth.n=512
smooth.roc.binormal <- function(roc, n) {
  df <- data.frame(sp=qnorm(roc$sp * ifelse(roc$percent, 1/100, 1)), se=qnorm(roc$se * ifelse(roc$percent, 1/100, 1)))
  df <- df[apply(df, 1, function(x) all(is.finite(x))),]
  if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
  model <- lm(sp~se, df)
  if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
  se <- qnorm(seq(0, 1, 1/(n-1)))
  sp <- predict(model, data.frame(se))

  return(list(sensitivities = pnorm(se) * ifelse(roc$percent, 100, 1),
              specificities = pnorm(sp) * ifelse(roc$percent, 100, 1),
              model = model))
}

# ???????ڣ???????ģ?ͣ???Ҫ???й????????ĵ?
se1 <- roc3[["FP"]][, 2]
sp1 <- roc3[["TP"]][, 2]
se1
sp1
# ???й????????ĵ?
se11 <- seq(se1[1], se1[2], length=10)
sp11 <- (sp1[1] - sp1[2])/(se1[1] - se1[2]) * (se11 - se1[2]) + sp1[2]
se12 <- seq(se1[2], se1[3], length=10)
sp12 <- (sp1[2] - sp1[3])/(se1[2] - se1[3]) * (se12 - se1[3]) + sp1[3]
se1 <- c(se11, se12)
sp1 <- c(sp11, sp12)
#
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t1 <- se2
sp_t1 <- sp2

# Ԥ??ģ??
se1 <- roc4[["FP"]][, 2]
sp1 <- roc4[["TP"]][, 2]
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t2 <- se2
sp_t2 <- sp2

# Ԥ??ģ??+????????
se1 <- roc5[["FP"]][, 2]
sp1 <- roc5[["TP"]][, 2]
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t3 <- se2
sp_t3 <- sp2

text1 <- paste0("Pathologic Stage model: AUC = ", sprintf("%.3f", roc3[["AUC"]][2]))
text2 <- paste0("Prognostic model: AUC = ", sprintf("%.3f", roc4[["AUC"]][2]))
text3 <- paste0("Combined model: AUC = ", sprintf("%.3f", roc5[["AUC"]][2]))

outfile <- "fig_roc_smooth_year3.pdf"
pdf(outfile, onefile=FALSE, width=15/2.54, height=15/2.54)
par(mar=c(4.5, 4.5, 2, 2))
cols <- c("red", "blue", "darkgreen")
ltys <- c(1, 1, 1, 1, 1)
lwds <- c(6, 6, 6, 6, 6)
legends <- c(text1, text2, text3)
plot(
    se_t1, sp_t1, type="l", col=cols[1], lty=ltys[1], lwd=lwds[1],
    xlim=c(0, 1), ylim=c(0, 1),
    xlab="1 - Specificity",
    ylab="Sensitivity"
)
polygon(c(se_t1, 1), c(sp_t1, 0), col="lightblue")
lines(se_t2, sp_t2, col=cols[2], lty=ltys[2], lwd=lwds[2])
lines(se_t3, sp_t3, col=cols[3], lty=ltys[3], lwd=lwds[3])
abline(0, 1, col="grey")
legend(
    0.9, 0.4, legend=legends, text.col=cols, adj=1, box.col=NA
)
text(x=0.6, y=0.42, labels="3 year AUC")
dev.off()

outfile <- "fig_roc_smooth_year3.tiff"
tiff(outfile, width=15, height=15, unit="cm", res=600, compression="lzw+p")
par(mar=c(4.5, 4.5, 2, 2))
cols <- c("red", "blue", "darkgreen")
ltys <- c(1, 1, 1, 1, 1)
lwds <- c(6, 6, 6, 6, 6)
legends <- c(text1, text2, text3)
plot(
    se_t1, sp_t1, type="l", col=cols[1], lty=ltys[1], lwd=lwds[1],
    xlim=c(0, 1), ylim=c(0, 1),
    xlab="1 - Specificity",
    ylab="Sensitivity"
)
polygon(c(se_t1, 1), c(sp_t1, 0), col="lightblue")
lines(se_t2, sp_t2, col=cols[2], lty=ltys[2], lwd=lwds[2])
lines(se_t3, sp_t3, col=cols[3], lty=ltys[3], lwd=lwds[3])
abline(0, 1, col="grey")
legend(
    0.9, 0.4, legend=legends, text.col=cols, adj=1, box.col=NA
)
text(x=0.6, y=0.42, labels="3 year AUC")
dev.off()




##################################################??5??ROC#######################################
# ??5????ROC
library(timeROC)
fit <- coxph(Surv(futime, fustat) ~ PathologicStage, data=bc)
# ????????
bc$pred <- predict(fit, newdata=bc)
roc3 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting='aalen',    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*5,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc3 <- roc3$AUC[2]
auc3


fit <- coxph(Surv(futime, fustat) ~ Risk, data=bc)
# Ԥ??ģ??
bc$pred <- predict(fit, newdata=bc)
roc4 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting='aalen',    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*5,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc4 <- roc4$AUC[2]
auc4


fit <- coxph(Surv(futime, fustat) ~ Risk + PathologicStage, data=bc)
# Ԥ??ģ??+????????
bc$pred <- predict(fit, newdata=bc)
roc5 <- timeROC(
    T=bc$futime,               #????ʱ??
    delta=bc$fustat,         #????????
    marker=bc$pred,     #Ԥ????��
    cause=1,                 #???Խ??ָ?ֵ??????????????
    weighting='aalen',    #Ȩ?ؼ??㷽????marginal??Ĭ??ֵ??????km????ɾʧ?ֲ?
    times=365*5,             #ʱ???㣬ѡȡ5??(60????)??8????????
    ROC=TRUE)
auc5 <- roc5$AUC[2]
auc5

#pdf(file="ROC-5??.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
plot(roc3$FP[, 2], roc3$TP[, 2], type="l", xlim=c(0,1), ylim=c(0,1),col=2, 
    xlab="False positive rate", ylab="True positive rate",
    lwd=6)
lines(roc4$FP[, 2], roc4$TP[, 2], type="l",col=3,lwd=6)
lines(roc5$FP[, 2], roc5$TP[, 2], type="l",col=4,lwd=6)
abline(0,1)
legend("bottomright", legend=c(
"Pathologic Stage model",
"Prognostic model",
"Combined model"),
lty=c(1, 1, 1, 1, 1), col=c(2, 3, 4), box.col="transparent",lwd=6)
# 1??ʵ?֣?2??????
text(x=0.6, y =0.2, labels ="5 year AUC")
#dev.off()
auc3 # ????????
auc4 # Ԥ??ģ??
auc5 # Ԥ??ģ??+????????


# ??Ϊƽ??????
# https://github.com/cran/pROC/blob/master/R/smooth.R
# smooth.method="binormal", smooth.n=512
smooth.roc.binormal <- function(roc, n) {
  df <- data.frame(sp=qnorm(roc$sp * ifelse(roc$percent, 1/100, 1)), se=qnorm(roc$se * ifelse(roc$percent, 1/100, 1)))
  df <- df[apply(df, 1, function(x) all(is.finite(x))),]
  if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
  model <- lm(sp~se, df)
  if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
  se <- qnorm(seq(0, 1, 1/(n-1)))
  sp <- predict(model, data.frame(se))

  return(list(sensitivities = pnorm(se) * ifelse(roc$percent, 100, 1),
              specificities = pnorm(sp) * ifelse(roc$percent, 100, 1),
              model = model))
}

# ???????ڣ???????ģ?ͣ???Ҫ???й????????ĵ?
se1 <- roc3[["FP"]][, 2]
sp1 <- roc3[["TP"]][, 2]
se1
sp1
# ???й????????ĵ?
se11 <- seq(se1[1], se1[2], length=10)
sp11 <- (sp1[1] - sp1[2])/(se1[1] - se1[2]) * (se11 - se1[2]) + sp1[2]
se12 <- seq(se1[2], se1[3], length=10)
sp12 <- (sp1[2] - sp1[3])/(se1[2] - se1[3]) * (se12 - se1[3]) + sp1[3]
se1 <- c(se11, se12)
sp1 <- c(sp11, sp12)
#
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t1 <- se2
sp_t1 <- sp2

# Ԥ??ģ??
se1 <- roc4[["FP"]][, 2]
sp1 <- roc4[["TP"]][, 2]
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t2 <- se2
sp_t2 <- sp2

# Ԥ??ģ??+????????
se1 <- roc5[["FP"]][, 2]
sp1 <- roc5[["TP"]][, 2]
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t3 <- se2
sp_t3 <- sp2

text1 <- paste0("Pathologic Stage model: AUC = ", sprintf("%.3f", roc3[["AUC"]][2]))
text2 <- paste0("Prognostic model: AUC = ", sprintf("%.3f", roc4[["AUC"]][2]))
text3 <- paste0("Combined model: AUC = ", sprintf("%.3f", roc5[["AUC"]][2]))

outfile <- "fig_roc_smooth_year5.pdf"
pdf(outfile, onefile=FALSE, width=15/2.54, height=15/2.54)
par(mar=c(4.5, 4.5, 2, 2))
cols <- c("red", "blue", "darkgreen")
ltys <- c(1, 1, 1, 1, 1)
lwds <- c(6, 6, 6, 6, 6)
legends <- c(text1, text2, text3)
plot(
    se_t1, sp_t1, type="l", col=cols[1], lty=ltys[1], lwd=lwds[1],
    xlim=c(0, 1), ylim=c(0, 1),
    xlab="1 - Specificity",
    ylab="Sensitivity"
)
polygon(c(se_t1, 1), c(sp_t1, 0), col="lightblue")
lines(se_t2, sp_t2, col=cols[2], lty=ltys[2], lwd=lwds[2])
lines(se_t3, sp_t3, col=cols[3], lty=ltys[3], lwd=lwds[3])
abline(0, 1, col="grey")
legend(
    0.9, 0.4, legend=legends, text.col=cols, adj=1, box.col=NA
)
text(x=0.6, y=0.42, labels="5 year AUC")
dev.off()

outfile <- "fig_roc_smooth_year5.tiff"
tiff(outfile, width=15, height=15, unit="cm", res=600, compression="lzw+p")
par(mar=c(4.5, 4.5, 2, 2))
cols <- c("red", "blue", "darkgreen")
ltys <- c(1, 1, 1, 1, 1)
lwds <- c(6, 6, 6, 6, 6)
legends <- c(text1, text2, text3)
plot(
    se_t1, sp_t1, type="l", col=cols[1], lty=ltys[1], lwd=lwds[1],
    xlim=c(0, 1), ylim=c(0, 1),
    xlab="1 - Specificity",
    ylab="Sensitivity"
)
polygon(c(se_t1, 1), c(sp_t1, 0), col="lightblue")
lines(se_t2, sp_t2, col=cols[2], lty=ltys[2], lwd=lwds[2])
lines(se_t3, sp_t3, col=cols[3], lty=ltys[3], lwd=lwds[3])
abline(0, 1, col="grey")
legend(
    0.9, 0.4, legend=legends, text.col=cols, adj=1, box.col=NA
)
text(x=0.6, y=0.42, labels="5 year AUC")
dev.off()


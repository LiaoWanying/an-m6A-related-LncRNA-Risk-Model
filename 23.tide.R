# Initialize settings
setwd("//Users/lwy/Desktop/重跑v2/data15/2")
options(stringsAsFactors=FALSE)
#wd <- "C:/data/expl/order/20210614766/"
#setwd(wd)


# ????????
infile <- "allRisk.txt"
clindf <- read.table(infile, header=TRUE, sep="\t", check.names=FALSE)
dim(clindf)
clindf[1:6,]
clindf1 <- clindf


infile <- "TIDE_output-jieduan.csv"
clindf <- read.csv(infile, header=TRUE, check.names=FALSE)
dim(clindf)
clindf[1:6,]
clindf2 <- clindf
colnames(clindf2)[1] <- "id"


# ȡ????
clindf <- merge(clindf1, clindf2, by="id", all=FALSE)
dim(clindf)
clindf[1:6,]

# ?鿴????????
sample_ids <- clindf[, 1]
sample_ids <- sample_ids[order(sample_ids)]
types <- substr(sample_ids, 14, 15)
table(types)
types <- ifelse(types >= 10, "normal", "tumor")
table(types)
# ɾ??????????
sample_ids <- sample_ids[types != "normal"]
length(sample_ids)
# ?鿴ͬһ???˶???????
patient_ids <- substr(sample_ids, 1, 12)
sample_ids[patient_ids %in% patient_ids[duplicated(patient_ids)]]
# ͬһ???˲?ͬ??????????????һ??
sample_ids <- sample_ids[!duplicated(patient_ids)]
length(sample_ids)


# ?ٲ?ͼ
library(ggplot2)
library(ggsci)
newdf1 <- clindf
newdf1 <- newdf1[newdf1[, "risk"] == "high",]
newdf1[, "Responder"] <- factor(
    newdf1[, "Responder"], levels=c("FALSE", "TRUE"),
    labels=c("Non-response", "Response")
)
table(newdf1[, "Responder"])
newdf1 <- newdf1[order(-newdf1[, "TIDE"]),]
newdf1[, "id2"] <- 1:nrow(newdf1)


# ??ɫ
colors <- c("blue", "red")

plot <- ggplot(newdf1, aes(x=id2, y=TIDE, fill=Responder)) +
    geom_bar(stat="identity", width=0.5) +
    scale_x_continuous("High", expand=c(0, 0)) +
    scale_y_continuous("TIDE prediction score") +
    scale_fill_manual(values=colors, guide=guide_legend(row=1, byrow=TRUE)) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y=element_line(),
        axis.text.y=element_text(),
        axis.title=element_text(),
        legend.title=element_blank(),
        legend.text=element_text(),
        legend.position=c(0.7, 0.8),
        panel.background=element_blank(),
        #panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
    )
plot

outfile <- "fig_barplot_waterfall_high.tiff"
tiff(outfile, width=20, height=10, units="cm", res=600, compression="lzw+p")
print(plot)
dev.off()

outfile <- "fig_barplot_waterfall_high.pdf"
pdf(outfile, width=20/2.54, height=10/2.54)
print(plot)
dev.off()


# Low
newdf1 <- clindf
newdf1 <- newdf1[newdf1[, "risk"] == "low",]
newdf1[, "Responder"] <- factor(
    newdf1[, "Responder"], levels=c("FALSE", "TRUE"),
    labels=c("Non-response", "Response")
)
table(newdf1[, "Responder"])
newdf1 <- newdf1[order(-newdf1[, "TIDE"]),]
newdf1[, "id2"] <- 1:nrow(newdf1)

colors <- c("blue", "red")

plot <- ggplot(newdf1, aes(x=id2, y=TIDE, fill=Responder)) +
    geom_bar(stat="identity", width=0.5) +
    scale_x_continuous("Low", expand=c(0, 0)) +
    scale_y_continuous("TIDE prediction score") +
    scale_fill_manual(values=colors, guide=guide_legend(nrow=1, byrow=TRUE)) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y=element_line(),
        axis.text.y=element_text(),
        axis.title=element_text(),
        legend.title=element_blank(),
        legend.text=element_text(),
        legend.position=c(0.7, 0.8),
        panel.background=element_blank(),
        #panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
    )
plot

outfile <- "fig_barplot_waterfall_low.tiff"
tiff(outfile, width=30, height=10, units="cm", res=600, compression="lzw+p")
plot
dev.off()

outfile <- "fig_barplot_waterfall_low.pdf"
pdf(outfile, width=30/2.54, height=10/2.54)
plot
dev.off()


# ?ٷֱ???ͼ
library(ggplot2)
newdf1 <- clindf
newdf1[, "Responder"] <- factor(
    newdf1[, "Responder"], levels=c("FALSE", "TRUE"),
    labels=c("Non-response", "Response")
)
tb <- table(newdf1[, "risk"], newdf1[, "Responder"])
tb
model <- fisher.test(tb)
model
pvalue <- model[["p.value"]]
pvalue <- ifelse(
    pvalue < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", pvalue))
)

newdf2 <- as.data.frame(table(newdf1[, "risk"], newdf1[, "Responder"]))
newdf3 <- as.data.frame(table(newdf1[, "risk"]))
colnames(newdf3)[2] <- "sum"
newdf4 <- merge(newdf2, newdf3, by="Var1")
newdf4[, "pct"] <- newdf4[, "Freq"] / newdf4[, "sum"]
get_labely <- function(indf, group, value, adjust=0.5) {
    newdf <- indf
    labely <- NULL
    groups <- levels(factor(newdf[, group]))
    for (i in 1:length(groups)) {
        subdf <- newdf[newdf[, group] == groups[i],]
        labely <- c(
            labely, 1 - (cumsum(subdf[, value]) - subdf[, value]*adjust)
        )
    }
    return(labely)
}
pos <- get_labely(newdf4, "Var1", "pct")
pcts <- paste0(round(newdf4[, "pct"] * 100, 1), "%")
colors <- c("blue", "red")

plot <- ggplot(newdf4, aes(x=Var1, y=pct, fill=Var2)) +
    geom_bar(stat="identity", width=0.7) +
    geom_text(aes(y=pos, label=pcts)) +
    scale_fill_manual(values=colors) +
    scale_x_discrete("") +
    scale_y_continuous("", expand=c(0, 0)) +
    ggtitle(paste0("Fisher's exact test, ", pvalue)) +
    theme_bw() +
    theme(
        legend.position="top",
        legend.title=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(color=NA, fill=NA),
        plot.title=element_text(hjust=0.5)
    )
plot

outfile <- "fig_barplot_percent.tiff"
tiff(outfile, width=10, height=10, unit="cm", res=600, compression="lzw+p")
plot
dev.off()

outfile <- "fig_barplot_percent.pdf"
pdf(outfile, width=10/2.54, height=10/2.54)
plot
dev.off()


# ??ͼ
library(ggplot2)
library(ggpubr)
newdf1 <- clindf
colors <- c("red", "blue")


plot <- ggplot(newdf1, aes(x=risk, y=TIDE, color=risk)) +
    # ??????????Сsize?????õ??Ĵ?Сoutlier.size?????????ߴ?ϸfatten
    geom_boxplot(
        aes(color=risk), notch=FALSE,
        outlier.size=0.6, outlier.shape=1, outlier.alpha=0,
        size=1.2, fatten=1
    ) +
    geom_jitter() +
    scale_color_manual(values=colors) +
    scale_x_discrete("") +
    scale_y_continuous("") +
    stat_compare_means(
        aes(
            label=ifelse(
                p < 0.001, "P < 0.001", paste0("P = ", ..p.format..)
            )
        ),
        # ????Pֵ????labex.x
        label.x=1.5, hjust=0.5
    ) +
    theme_bw() +
    theme(
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=0, hjust=0.5),
        axis.text.y=element_text(),
        axis.line=element_line(),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA)
    )
plot

outfile <- "fig_boxplot.tiff"
tiff(outfile, width=10, height=10, unit="cm", res=600, compression="lzw+p")
plot
dev.off()

outfile <- "fig_boxplot.pdf"
pdf(outfile, width=10/2.54, height=10/2.54)
plot
dev.off()


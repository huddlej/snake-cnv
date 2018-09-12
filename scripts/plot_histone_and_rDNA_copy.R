library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.file <- args[2]

df <- read.table(input.file, header=TRUE)
wide.df <- dcast(df, sample ~ region, value.var="copy_number")

pdf(output.file, width=8, height=6)
ggplot(wide.df, aes(x=`chr2L:21420129-21420657`, y=`chrUn_CP007120v1:43484-45478`)) + geom_point() + theme_classic() + labs(x="Histone copy number", y="rDNA copy number") + xlim(c(0, 600)) + ylim(c(0, 600)) + geom_smooth(method="lm")
dev.off()

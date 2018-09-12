library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.file <- args[2]

data <- read.table(input.file, header=TRUE)

pdf(output.file, width=8, height=6)
ggplot(data, aes(x=type, y=mean_cn)) + geom_violin() + geom_boxplot(width=0.2) + theme_classic() + labs(x="CNV type", y="Copy number") + ylim(c(0, 10)) + geom_hline(yintercept=2) + geom_hline(yintercept=1, colour="red") + geom_hline(yintercept=3, colour="red")
dev.off()

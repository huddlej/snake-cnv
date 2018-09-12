library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.file <- args[2]

data <- read.table(input.file, header=TRUE)

pdf(output.file, width=8, height=6)
ggplot(data, aes(x=diploid_status, y=copy_number)) + geom_violin() + geom_boxplot(width=0.2) + theme_classic() + labs(x="Copy number status", y="Copy number") + ylim(c(0, 10)) + geom_hline(yintercept=2)
dev.off()

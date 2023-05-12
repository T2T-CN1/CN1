args <- commandArgs(trailingOnly = TRUE)
if ( is.null(args[1])) {
	stop("Need sample.histo !")
}
print(args)

output <- paste(sep = '',args[1],".pdf")
kmer <- read.table(args[1])
kmer <- subset(kmer, V1 >=1 & V1 <=500) #对频数范围5-500的数据进行绘制 
Frequency <- kmer$V1
Number <- kmer$V2
pdf(output)
plot(Frequency, Number, type = 'l', col = 'blue')
dev.off()

library(dplyr)

Args <- commandArgs(TRUE)

all_reads <- read.table(Args[1], header = F, sep = "\t")
exon_reads <- read.table(Args[2], header = F, sep = "\t")

filter_reads <- all_reads[which(! all_reads$V4 %in% exon_reads$V4), ]
filter_reads_2 <- filter_reads %>%
    distinct(V4, .keep_all = T)

## over-gene reads
write.table(filter_reads_2, Args[3], sep = "\t",
    row.names = F, col.names = F, quote = F)

## intergenic reads except over-gene reads
itg <- read.table(Args[4], header = F, sep = "\t")
itg_reads <- itg[which(! itg$V4 %in% filter_reads_2$V4), ]
write.table(itg_reads, Args[5], sep = "\t",
    row.names = F, col.names = F, quote = F)
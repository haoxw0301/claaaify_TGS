Args <- commandArgs(TRUE) # nolint

fusion_reads <- read.table(Args[1], header = F, sep = "\t")
wi <- read.table(Args[2], header = F, sep = "\t")
wo <- read.table(Args[3], header = F, sep = "\t")

read_through <- fusion_reads[
    which(! fusion_reads$V4 %in% rbind(wi, wo)$V4), ]

write.table(read_through, Args[4], sep = "\t",
    row.names = F, col.names = F, quote = F)
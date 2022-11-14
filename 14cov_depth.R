library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
####---- load data -----#####
Args <- commandArgs(TRUE) # nolint
cov_list <- read.table(Args[1], header = F, sep = "\t")
reads_list <- read.table(Args[2], header = F, sep = "\t")

all_reads <- read.table(paste(Args[3], ".bed12", sep = ""),
                header = F, sep = "\t")

### ---- calculate cov_depth ----###
cov_depth <- function(df) {
    cov_length <- sum(df$V5) # nolint
    cov_ratio <- cov_length / sum(df$V6) # nolint

    df$depth <- as.numeric(df$V4) * as.numeric(df$V5) / df$V6 # nolint
    mean_depth <- mean(df$depth) # nolint
    return(list(cov_length = cov_length, cov_ratio = cov_ratio,
                mean_depth = mean_depth))
}

filelist <- c()
cov_length_list <- c()
cov_ratio_list <- c()
mean_depth_list <- c()
for (i in 1:nrow(cov_list)) { # nolint
    file <- cov_list[i, 1]
    df <- read.table(file, header = F, sep = "\t")
    df$depth <- as.numeric(df$V4) * as.numeric(df$V5) / df$V6 # nolint
    filelist[i] <- file
    cov_length_list[i] <- cov_depth(df)$cov_length
    cov_ratio_list[i] <- cov_depth(df)$cov_ratio
    mean_depth_list[i] <- cov_depth(df)$mean_depth
}
### calculate the sum length of the exon of each read
rna_length <- function(row) {
    exon_length <- strsplit(row[11], ",")[[1]]

    RNA <- sum(as.numeric(exon_length)) # nolint
    return(RNA)
}

### calculate the statistics of each group of reads
reads_statistics <- function(df) {
    df$rna_length <- apply(df, 1, rna_length)
    df$cover_length <- as.numeric(df$V3) - as.numeric(df$V2) # nolint
    names(df)[10] <- "exon_num"
    mean_exon_length <- sum(df$rna_length) / sum(as.numeric(df$exon_num)) # nolint
    mean_cover_length <- sum(df$cover_length) / nrow(df) # nolint
    mean_RNA_length <- sum(df$rna_length) / nrow(df) # nolint

    return(list(mean_exon_length = mean_exon_length,
                mean_cover_length = mean_cover_length,
                mean_RNA_length = mean_RNA_length))
}

#### res
cov_res <- data.frame(
    file = filelist,
    cov_length = cov_length_list,
    cov_ratio = cov_ratio_list,
    mean_depth = mean_depth_list
)

cov_res <- separate(cov_res, col = file,
                into = c("file_1", "file_2"),
                sep = "/") %>%
            separate(col = file_2,
                into = c("file_1_1", "file_1_2"),
                sep = "_cov") %>%
            separate(col = file_1_1,
                into = c("file_1_1_1", "file_1_1_2"),
                sep = paste(as.character(Args[3]), "_", sep = ""))
cov_res_2 <- cov_res[, c(3, 5:7)]

#### depth
mel <- c() ## mean exon length
mcl <- c() ## mean cover length
mrl <- c() ## mean RNA length
rc_ratio <- c() ## reads count ratio

df_res <- data.frame(
    V1 = "V1", V2 = 0, V3 = 0, V4 = "V4", V5 = 0, V6 = "-",
    V7 = 0, V8 = 0, V9 = "V9", V10 = 0, V11 = "V11", V12 = "V12", file = "type"
)
#########------ calculate the length distribution and mean length ------######
for(i in 1:nrow(reads_list)) { # nolint
    file <- reads_list[i, 1]
    df <- read.table(file, header = F, sep = "\t")
    df$file <- file
    print(file)
    # print(nrow(df))
    df_res <- rbind(df_res, df)
    head(df_res)
    mel[i] <- reads_statistics(df)$mean_exon_length
    mcl[i] <- reads_statistics(df)$mean_cover_length
    mrl[i] <- reads_statistics(df)$mean_RNA_length
    rc_ratio[i] <- nrow(df) / nrow(all_reads)
}

df_res <- df_res[-1, ]

#### length distribution ------
df_res$rna_length <- apply(df_res, 1, rna_length)
df_res$cover_length <- as.numeric(df_res$V3) - as.numeric(df_res$V2)

df_res <- separate(df_res, col = file,
                into = c("file_1", "file_2"),
                sep = "/") %>%
            separate(col = file_2,
                into = c("file_1_1", "file_1_2"),
                sep = "\\.") %>%
            separate(col = file_1_1,
                into = c("file_1_1_1", "file_1_1_2"),
                sep = paste(as.character(Args[3]), "_", sep = ""))
df_res_2 <- df_res[, c(1:12, 15, 17, 18)]
names(df_res_2)[13] <- "type"

df_res_2$type <- factor(df_res_2$type,
        levels = c("genic",
            "coding",
            "lncRNA",
            "well-spliced",
            "ill-spliced",
            "intron-alone",
            "beyond-gene",
            "beyond-gene_5E",
            "beyond-gene_3E",
            # "beyond-gene_BE",
            "intergenic",
            "over-gene",
            "proximal_intergenic",
            "distal_intergenic"),
            ordered = TRUE)


df_res_3 <- melt(df_res_2,
            id.vars = colnames(df_res_2)[1:13])

pdf(paste("classify_res_2/", Args[3], "_cov_rna_length.pdf", sep = ""),
    width = 15, height = 5)
ggplot(data = df_res_3, aes(x = type, y = value, fill = variable)) +
    #geom_violin(color = NA) +
    geom_boxplot(width = 0.5, outlier.colour = NA) +
    scale_y_log10() +
    theme_bw() +
    scale_fill_manual(values = c("green", "orange")) +
    labs(x = "", y = "length", fill = "") +
    theme(axis.title = element_text(size = 20),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(hjust = 1, angle = 60, size = 16),
    legend.text = element_text(size = 14))
dev.off()
### splicing times ----

df_res_2$splicing_times <- as.numeric(df_res_2$V10) - 1

itg_res <- df_res_2[
    which(df_res_2$type %in% c(
            "over-gene",
            "proximal_intergenic",
            "distal_intergenic")), ]

pdf(paste("classify_res_2/", Args[3], "_cov_splicing_times.pdf", sep = ""),
    width = 10, height = 5)
ggplot(data = itg_res, aes(x = splicing_times, fill = type)) +
    geom_bar(position = "dodge") +
    theme_bw() +
    scale_y_log10() +
    scale_fill_manual(values = c("#3d435d", "royalblue4", "#1b2458")) +
    labs(x = "splicing times", y = "count", fill = "") +
    theme(axis.title = element_text(size = 20),
    axis.text = element_text(size = 14, color = "black"),
    # axis.text.x = element_text(hjust = 1, angle = 60, size = 16),
    legend.text = element_text(size = 14))
dev.off()

print(
    paste(
        "itg splicing ratio = ",
    nrow(itg_res[which(itg_res$splicing_times != 0), ]) / nrow(itg_res))
)
#### mean length of all reads -----
res <- data.frame(
    file = reads_list$V1,
    mean_exon_length = mel,
    mean_cover_length = mcl,
    mean_RNA_length = mrl,
    reads_count_ratio = rc_ratio
)

res <- separate(res, col = file,
                into = c("file_1", "file_2"),
                sep = "/") %>%
            separate(col = file_2,
                into = c("file_1_1", "file_1_2"),
                sep = "\\.") %>%
            separate(col = file_1_1,
                into = c("file_1_1_1", "file_1_1_2"),
                sep = paste(as.character(Args[3]), "_", sep = ""))

res_2 <- res[, c(3, 5:8)]

res_3 <- merge(res_2, cov_res_2, by = "file_1_1_2")
res_3$type <- as.character(Args[3])

res_3$reads_count_ratio <- round(res_3$reads_count_ratio, 4)
res_3$cov_length <- round(res_3$cov_length, 4)
res_3$mean_depth <- round(res_3$mean_depth, 4)

write.table(res_3, Args[4],
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
####---- level a ----
# cov_a <- res_3[which(res_3$file_1_1_2 %in%
#                     c("genic", "intergenic", "fusion")), c(1, 7)]

cov_a <- rbind(
    res_3[which(res_3$file_1_1_2 == "genic"), c(1, 7)],
    res_3[which(res_3$file_1_1_2 == "intergenic"), c(1, 7)]
)
cov_a <- rbind(cov_a,
            data.frame(file_1_1_2 = "other",
            cov_ratio = 1 - sum(cov_a$cov_ratio)))

cols_a <- c("genic" = "red", "intergenic" = "blue",
            "no_cover" = "grey")

cov_a$piepercent <- round(100 * cov_a$cov_ratio, 1)

cov_a$label <- paste(
    cov_a$file_1_1_2, " ", cov_a$piepercent, "%", sep = ""
)

cov_a$file_1_1_2 <- factor(cov_a$file_1_1_2,
        levels = c("genic", "intergenic", "no-cover"),
        ordered = TRUE)

pdf(paste("classify_res_2/", Args[3], "_cov_ratio.pdf", sep = ""),
    width = 10, height = 5)
pie(cov_a$cov_ratio, labels = NA, col = cols_a, cex = 1,
    main = "cover ratio of all reads")
legend("right", legend = cov_a$label, fill = cols_a, cex = 1)
dev.off()


print("04")
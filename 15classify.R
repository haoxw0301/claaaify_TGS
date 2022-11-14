####---- required packages ----
library(tidyr)
library(dplyr)
library(ggplot2)

Args <- commandArgs(TRUE) # nolint

####---- annotation data ----####
reads_num <- read.table(
    paste("classify_res_2/", Args[1], "_reads_num", sep = ""))

reads_num_2 <- separate(reads_num,
                        col = V2,
                        into = c("dic", "file"),
                        sep = paste(Args[1], "_", sep = "")) %>%
                        separate(col = file,
                        into = c("type", "bed12"),
                        sep = "\\.")
reads_num_2 <- reads_num_2[, c(1, 3)]

####---- level a ----
reads_num_a <- reads_num_2[
    which(reads_num_2$type %in% c("genic", "intergenic")), ]

cols_a <- data.frame(
    type = c("genic", "intergenic"),
    color = c("red", "blue"))

reads_num_a <- merge(reads_num_a, cols_a, by = "type")

reads_num_a$piepercent <- round(100 * reads_num_a$V1 / sum(reads_num_a$V1), 1)

reads_num_a$label <- paste(
    reads_num_a$type, " ", reads_num_a$piepercent, "%", sep = "")

pdf(paste("classify_res_2/", Args[1], "_reads_num_a.pdf", sep = ""),
    width = 5, height = 7)
pie(reads_num_a$V1, labels = NA, col = reads_num_a$color, cex = 1,
    main = "All reads")
legend("bottom", legend = reads_num_a$label, fill = reads_num_a$color, cex = 1)
dev.off()

####---- level b1, genic reads----
reads_num_b1 <- reads_num_2[
    which(reads_num_2$type %in%
        c("coding", "lncRNA")), ]
cols_b1 <- data.frame(
    type = c("coding", "lncRNA"),
    color = c("brown4", "brown2"))
reads_num_b1 <- merge(reads_num_b1, cols_b1, by = "type")
reads_num_b1$piepercent <- round(
    100 * reads_num_b1$V1 / sum(reads_num_a$V1), 2)
reads_num_b1$label <- paste(
    reads_num_b1$type, " ", reads_num_b1$piepercent, "%", sep = "")

reads_num_b1 <- reads_num_b1[which(reads_num_b1$piepercent > 0), ]
pdf(paste("classify_res_2/", Args[1], "_reads_num_b1.pdf", sep = ""),
    width = 5, height = 7)
pie(reads_num_b1$V1,
    labels = reads_num_b1$type,
    col = reads_num_b1$color,
    main = "genic reads", cex = 1)
legend("bottom", legend = reads_num_b1$label, fill = reads_num_b1$color, cex = 1)
dev.off()

####---- level b2, intergenic reads ----
reads_num_b2 <- reads_num_2[
    which(reads_num_2$type %in%
        c("proximal_intergenic", "distal_intergenic", "over-gene")), ]

cols_b2 <- data.frame(
    type = c("proximal_intergenic", "distal_intergenic", "over-gene"),
    color = c("royalblue4", "#1b2458", "#3d435d"))

reads_num_b2 <- merge(reads_num_b2, cols_b2, by = "type")

reads_num_b2$piepercent <- round(
    100 * reads_num_b2$V1 / sum(reads_num_a$V1), 2)

reads_num_b2$label <- paste(
    reads_num_b2$type, " ", reads_num_b2$piepercent, "%", sep = "")

reads_num_b2 <- reads_num_b2[which(reads_num_b2$piepercent > 0), ]
pdf(paste("classify_res_2/", Args[1], "_reads_num_b2.pdf", sep = ""),
    width = 5, height = 7)
pie(reads_num_b2$V1,
    labels = reads_num_b2$type,
    col = reads_num_b2$color,
    main = "intergenic reads", cex = 1)
legend("bottom", legend = reads_num_b2$label, fill =reads_num_b2$color, cex = 1)
dev.off()


####---- level c1, coding reads ----
reads_num_c1 <- reads_num_2[
    which(reads_num_2$type %in%
        c("well-spliced",
        "ill-spliced",
        "intron-alone",
        "beyond-gene")), ]

cols_c1 <- data.frame(
    type = c("well-spliced", "ill-spliced", "intron-alone", "beyond-gene"),
    color = c("violetred4", "violetred2", "violetred1", "violetred3"))

reads_num_c1 <- merge(reads_num_c1, cols_c1, by = "type")

reads_num_c1$piepercent <- round(
    100 * reads_num_c1$V1 / sum(reads_num_a$V1), 2)

reads_num_c1$label <- paste(
    reads_num_c1$type, " ", reads_num_c1$piepercent, "%", sep = "")

# reads_num_c1 <- reads_num_c1[which(reads_num_c1$piepercent > 0), ]

pdf(paste("classify_res_2/", Args[1], "_reads_num_c1.pdf", sep = ""),
    width = 5, height = 7)
pie(reads_num_c1$V1,
    labels = reads_num_c1$type,
    col = reads_num_c1$color,
    main = "coding reads", cex = 1)
legend("bottom", legend = reads_num_c1$label, fill = reads_num_c1$color, cex = 1)
dev.off()

####---- level d1, chimeric transcripts ----
reads_num_d1 <- reads_num_2[
    which(reads_num_2$type %in%
        c("beyond-gene_5E",
            "beyond-gene_3E",
            "beyond-gene_BE")), ]

cols_d1 <- data.frame(
    type = c("beyond-gene_5E", "beyond-gene_3E", "beyond-gene_BE"),
    color = c("darkorange4", "darkorange2", "darkorange1"))

reads_num_d1 <- merge(reads_num_d1, cols_d1, by = "type")

reads_num_d1$piepercent <- round(
    100 * reads_num_d1$V1 / sum(reads_num_a$V1), 3)

reads_num_d1$label <- paste(
    reads_num_d1$type, " ", reads_num_d1$piepercent, "%", sep = "")

reads_num_d1 <- reads_num_d1[which(reads_num_d1$piepercent > 0), ]

pdf(paste("classify_res_2/", Args[1], "_reads_num_d1.pdf", sep = ""),
    width = 5, height = 7)
pie(reads_num_d1$V1,
    labels = reads_num_d1$type,
    col =  reads_num_d1$color,
    main = "beyond-gene transcripts", cex = 1)
legend("bottom", legend = reads_num_d1$label, fill = reads_num_d1$color, cex = 1)
dev.off()
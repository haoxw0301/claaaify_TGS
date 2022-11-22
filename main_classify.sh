#!/bin/bash
i=${1}
i=$(basename ${i} ".bed12")  
####====== genic reads ======####
## protein coding genes, lncRNA 
    cat  ~/haoxw/reference/hg38/v42/protein_coding/transcript.gtf ~/haoxw/reference/hg38/v42/lncRNA/transcript.gtf | bedtools intersect -a ${i}.bed12 -b - -f 0.5 -wa | uniq > classify_tmp/${i}_genic.bed12
    
    bedtools bed12tobed6 -i classify_tmp/${i}_genic.bed12 > classify_tmp/${i}_genic.bed6

    cat  ~/haoxw/reference/hg38/v42/protein_coding/transcript.gtf ~/haoxw/reference/hg38/v42/lncRNA/transcript.gtf |bedtools intersect -a classify_tmp/${i}_genic.bed6 -b - | uniq > classify_tmp/${i}_genic_transcript.bed6

    Rscript 01filter_reads.R classify_tmp/${i}_genic.bed12 classify_tmp/${i}_genic_transcript.bed6 classify_res_2/${i}_genic.bed12

    ls  classify_res_2/${i}_genic.bed12 > classify_res_2/${i}_reads_list

    ## coverage
    bedtools coverage  -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_genic.bed12 > classify_res_2/${i}_genic_cov.txt
####----- protein coding reads -----####
    bedtools intersect -a ${i}.bed12 -b ~/haoxw/reference/hg38/v42/protein_coding/transcript_merge.bed -f 0.5  -wa | uniq > classify_tmp/${i}_coding.bed12
    
    bedtools bed12tobed6 -i classify_tmp/${i}_coding.bed12 > classify_tmp/${i}_coding.bed6
  
    bedtools intersect -a classify_tmp/${i}_coding.bed6 -b ~/haoxw/reference/hg38/v42/protein_coding/transcript_merge.bed -wa | uniq > classify_tmp/${i}_coding_transcript.bed6 ## isoform overlap with gene
  
    Rscript 01filter_reads.R classify_tmp/${i}_coding.bed12 classify_tmp/${i}_coding_transcript.bed6 classify_res_2/${i}_coding.bed12
    ### coverage
    bedtools coverage  -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_coding.bed12 > classify_res_2/${i}_coding_cov.txt

    ls classify_res_2/${i}_coding.bed12 >> classify_res_2/${i}_reads_list
####-----intron alone reads -----
    bedtools intersect -a ${i}.bed12 -b ~/haoxw/reference/hg38/v42/protein_coding/intron.bed -f 0.99  -wa | uniq > classify_res_2/${i}_intron-alone.bed12 ## reads fully mapped to a single intron

    ### coverage
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_intron-alone.bed12 > classify_res_2/${i}_intron-alone_cov.txt

    ls classify_res_2/${i}_intron-alone.bed12 >> classify_res_2/${i}_reads_list
####---- well and ill spilced reads ----
### reads mapped to both introns and exons
### reads mapped to intronless genes are classified as well-spliced reads
    ## intronless gene reads
    bedtools intersect -a classify_res_2/${i}_coding.bed12 -b ~/haoxw/reference/hg38/v42/protein_coding/transcript_intronless.bed -f 0.9 -wa | uniq > classify_res_2/${i}_intronless.bed12

    ## other genes
    
    bedtools bed12tobed6 -i classify_res_2/${i}_coding.bed12 > classify_tmp/${i}_coding_exon.bed6
    bedtools intersect -a classify_tmp/${i}_coding_exon.bed6 -b ~/haoxw/reference/hg38/v42/protein_coding/intron.bed -wo -s | uniq > classify_tmp/${i}_intersect_intron.txt ## 1
    ## exon reads intersect with exon
    bedtools intersect -a classify_tmp/${i}_coding_exon.bed6 -b ~/haoxw/reference/hg38/v42/*/exon_merge.bed -wo | uniq > classify_tmp/${i}_intersect_exon.txt ## 2

    ### intron and exon reads
    bedtools intersect -a ${i}.bed12 -b ~/haoxw/reference/hg38/v42/protein_coding/intron.bed -s -wa | uniq > classify_tmp/${i}_intersect_intron.bed12 ## 3
    bedtools intersect -a ${i}.bed12 -b ~/haoxw/reference/hg38/v42/protein_coding/exon_merge.bed -s -wa | uniq > classify_tmp/${i}_intersect_exon.bed12 ## 4
    
    Rscript 03splice_class.R classify_tmp/${i}_intersect_intron.txt classify_tmp/${i}_intersect_exon.txt classify_tmp/${i}_intersect_intron.bed12  classify_tmp/${i}_intersect_exon.bed12 classify_res_2/${i}_well-spliced.bed12 classify_res_2/${i}_poor-spliced.bed12 classify_res_2/${i}_non-spliced.bed12 classify_res_2/${i}_length_density.pdf

    cat classify_res_2/${i}_intronless.bed12  >> classify_res_2/${i}_well-spliced.bed12
    ## coverage
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_well-spliced.bed12 > classify_res_2/${i}_well-spliced_cov.txt

    ls classify_res_2/${i}_well-spliced.bed12 >> classify_res_2/${i}_reads_list

####------ ill spliced reads ------
## short exon reads
    bedtools intersect -a ${i}.bed12 -b ~/haoxw/reference/hg38/v42/2_exon_merge.bed -f 0.99  -wa | uniq > classify_tmp/${i}_exon.bed12 ## reads fully mapped to a single exon
    
    cat classify_res_2/${i}_poor-spliced.bed12 classify_res_2/${i}_non-spliced.bed12 classify_tmp/${i}_exon.bed12 > classify_res_2/${i}_ill-spliced.bed12

    ### coverage
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_ill-spliced.bed12 > classify_res_2/${i}_ill-spliced_cov.txt

    ls classify_res_2/${i}_ill-spliced.bed12 >> classify_res_2/${i}_reads_list

####---- Beyond gene reads ------#### 

    bedtools intersect -a classify_res_2/${i}_coding.bed12 -b ~/haoxw/reference/hg38/v42/protein_coding/non_gene.bed -f 0.01 -wa | uniq > classify_res_2/${i}_beyond-gene.bed12

    bedtools intersect -a classify_res_2/${i}_beyond-gene.bed12 -b ~/haoxw/reference/hg38/v42/protein_coding/transcript_merge.bed -wo  | uniq > classify_res_2/${i}_beyond-gene.txt

    Rscript 16beyond_gene.R classify_res_2/${i}_beyond-gene.txt classify_res_2/${i}_beyond-gene_5E.bed12 classify_res_2/${i}_beyond-gene_3E.bed12 
    
    ### coverage
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_beyond-gene.bed12 > classify_res_2/${i}_beyond-gene_cov.txt
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_beyond-gene_5E.bed12 > classify_res_2/${i}_beyond-gene_5E_cov.txt
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_beyond-gene_3E.bed12 > classify_res_2/${i}_beyond-gene_3E_cov.txt
    # bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_beyond-gene_BE.bed12 > classify_res_2/${i}_beyond-gene_BE_cov.txt

    ls classify_res_2/${i}_beyond-gene.bed12 >> classify_res_2/${i}_reads_list
    ls classify_res_2/${i}_beyond-gene_5E.bed12 >> classify_res_2/${i}_reads_list
    ls classify_res_2/${i}_beyond-gene_3E.bed12 >> classify_res_2/${i}_reads_list
    # ls classify_res_2/${i}_beyond-gene_BE.bed12 >> classify_res_2/${i}_reads_list
############----- lncRNA reads -----#################################

    bedtools intersect  -a ${i}.bed12 -b  ~/haoxw/reference/hg38/v42/lncRNA/transcript.gtf -f 0.5 -wa | uniq > classify_tmp/${i}_lncRNA.bed12
    
    bedtools subtract -A -a classify_tmp/${i}_lncRNA.bed12 -b ~/haoxw/reference/hg38/v42/protein_coding/intron.bed | uniq > classify_res_2/${i}_lncRNA.bed12

    ### coverage
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_lncRNA.bed12 > classify_res_2/${i}_lncRNA_cov.txt

    ls classify_res_2/${i}_lncRNA.bed12 >> classify_res_2/${i}_reads_list

#####################################################################
####====== intergenic reads ====####
#### reads have no overlap with above transcripts
  
    bedtools intersect -a ${i}.bed12 -b ~/haoxw/reference/hg38/v42/intergenic.bed -f 0.5 -wa | uniq > classify_res_2/${i}_intergenic.bed12

    ## over-gene reads
    cat ~/haoxw/reference/hg38/v42/protein_coding/transcript.gtf ~/haoxw/reference/hg38/v42/lncRNA/transcript.gtf | bedtools intersect -a ${i}.bed12 -b - -F 1 -wa | uniq > classify_tmp/${i}_over-gene.bed12

    bedtools bed12tobed6 -i classify_tmp/${i}_over-gene.bed12 > classify_tmp/${i}_over-gene.bed6
    cat ~/haoxw/reference/hg38/v42/protein_coding/transcript.gtf ~/haoxw/reference/hg38/v42/lncRNA/transcript.gtf | bedtools intersect -a classify_tmp/${i}_over-gene.bed6 -b - -wa | uniq > classify_tmp/${i}_over-gene_txp.bed6


    # R
    Rscript 13over_gene.R classify_tmp/${i}_over-gene.bed12 classify_tmp/${i}_over-gene_txp.bed6 classify_res_2/${i}_over-gene.bed12 classify_res_2/${i}_intergenic.bed12 classify_tmp/${i}_intergenic.bed12

    bedtools intersect -a classify_tmp/${i}_intergenic.bed12 -b ~/haoxw/reference/hg38/v42/distal_intergenic.bed -f 0.5 -wa | uniq > classify_res_2/${i}_distal_intergenic.bed12

    bedtools intersect -a classify_tmp/${i}_intergenic.bed12 -b ~/haoxw/reference/hg38/v42/proximal_intergenic.bed -f 0.5 -wa | uniq  > classify_res_2/${i}_proximal_intergenic.bed12

    rm classify_res_2/${i}_intergenic.bed12
    cat classify_res_2/${i}_over-gene.bed12 classify_res_2/${i}_distal_intergenic.bed12 classify_res_2/${i}_proximal_intergenic.bed12 > classify_res_2/${i}_intergenic.bed12
    ## coverage
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_intergenic.bed12 > classify_res_2/${i}_intergenic_cov.txt
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_over-gene.bed12 > classify_res_2/${i}_over-gene_cov.txt
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_proximal_intergenic.bed12 > classify_res_2/${i}_proximal_intergenic_cov.txt
    bedtools coverage -a ~/haoxw/reference/hg38/hg38_windows_100k.bed -b classify_res_2/${i}_distal_intergenic.bed12 > classify_res_2/${i}_distal_intergenic_cov.txt


    ls classify_res_2/${i}_intergenic.bed12 >> classify_res_2/${i}_reads_list
    ls classify_res_2/${i}_over-gene.bed12 >> classify_res_2/${i}_reads_list
    ls classify_res_2/${i}_proximal_intergenic.bed12 >> classify_res_2/${i}_reads_list
    ls classify_res_2/${i}_distal_intergenic.bed12 >> classify_res_2/${i}_reads_list


#####################################################################
#########===== coverage and depth analysis ========##################
    ls  classify_res_2/${i}*_cov.txt > classify_res_2/${i}_cov.list

    Rscript 14cov_depth.R classify_res_2/${i}_cov.list classify_res_2/${i}_reads_list ${i} classify_res_2/${i}_cov_depth_total.txt 
#####################################################################
rm classify_res_2/${i}_reads_num

for file in `cat classify_res_2/${i}_reads_list`
do
    wc -l $file >> classify_res_2/${i}_reads_num
done

Rscript 15classify.R ${i}   

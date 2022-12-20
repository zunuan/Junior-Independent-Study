library(biomaRt)
listEnsembl()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

chrom=c(14,7)
HUGO_list <- getBM(attributes= "hgnc_symbol", filters=c("chromosome_name"),values=list(chrom), mart=ensembl)

TRAV_list <- HUGO_list[grepl('TRAV', col$hgnc_symbol), ]
TRAJ_list <- HUGO_list[grepl('TRAJ', col$hgnc_symbol), ]

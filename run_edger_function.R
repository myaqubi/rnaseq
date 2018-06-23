library(edgeR)
library(biomaRt)
library(readr)
library(tidyverse)
library(magrittr)


run_edger <- function(counts_matrix, genome_name, groups_info, sample_threshold = 1, file_name){
  message("Showing first few rows of counts input")
  print(head(counts_matrix)) # for checking
  
  stopifnot(!is.null(row.names(counts_matrix))) # counts input needs to have rownames
  if(!grepl("EN", row.names(counts_matrix)[1])) message("WARNING : Rownames do not seem to be ensembl IDs") # rownames are not ensembl id
  counts = as.matrix(counts_matrix)
  
  makemart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset=genome_name   , host="www.ensembl.org")
  group = groups_info
  if(ncol(counts_matrix)!=length(group)) message("ERROR : Check the counts input") # counts input or groups info is not correct
  stopifnot(ncol(counts_matrix)==length(group))
  
  dim(counts)
  cds <- DGEList( counts , group = group )
  cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= sample_threshold, ]
  cds <- calcNormFactors( cds )
  cds <- estimateCommonDisp( cds )
  cds <- estimateTagwiseDisp( cds )
  etest <- exactTest(cds)
  top <- topTags(etest, n=nrow(cds$counts))$table
  de <- rownames(top)
  annotations <- getBM(attributes=c('ensembl_gene_id', 'description', 'external_gene_name'), filters = 'ensembl_gene_id', values = de, mart = makemart)
  annotations <- annotations[match(de, annotations$ensembl_gene_id), ]
  values <- as.data.frame(head(top, n=length(de)))
  dedata <- cbind(values, annotations)
  write_tsv(dedata, file_name)
  return (dedata)
  
}

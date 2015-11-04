#!/usr/bin/R
#Script to analyse Data from project of PGB

##Change these variables accordingly
cwd <- "/home/lluis/Documents/master/PGB/exercises/"
tissue <-  "Cerebellum" # Is used to search for the tissue
# if an error occurs explore the names(raw_data) and
# store the tissue_n the number of column of your tissue
raw_file <- "GTEx_Analysis_RNA-seq_selected_tissues.txt"

# Read the data and finds the tissue column
setwd(cwd)
col_classes <- c("character", "character", rep("numeric", 19))
raw_data <- read.table(raw_file, header=TRUE, sep="\t", nrows = 5278L, skip = 1)
tissue_n <-  grep(tissue, names(raw_data), ignore.case=TRUE)
tissue_name <- names(raw_data)[tissue_n]

# To connect to Biomart to retrieve data we can use these library 
library(biomaRt) 
# To install you need to do 
#source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

## Define functions needed for what we want to do
select_gene <- function(row){
  # Keep genest that have at least 1 tissue with 1 or more rpkm
  val = sum(row[-c(1,2)]>1)
  if (val>=1){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

type_gene <- function(row, tissue){
# Return the tag of the gene
  val_all = sum(row[-c(1,2)]>1)
  val_tissue = row[tissue] > 1
  if (val_all==19){
    return("Housekeeping")
  } else if (val_all<=3 & val_tissue){
    return("Tissue-specific")
  } else {
    return("Not-special")
  }
}

get_genes <- function(gene){
  # Get the Ensembl identifier for each gene by spliting the version
  g <- strsplit(as.character(gene), split =".", fixed=TRUE)
  return(g[[1]][1])
}


# Select which genes are above threshold
keep <- apply(raw_data, 1, select_gene)
data <- raw_data[keep,]

# Remove the version of the gene
genes_id <- sapply(data$Gene, get_genes)
data$Gene <- genes_id

# Classify the remaining genes
tag <-  as.factor(apply(data, 1, type_gene, tissue=tissue_n))
dat <- cbind(data, tag)

# Select the tissue genes
col <- c("Gene", "Protein", tissue_name, "tag")
data_c <- dat[, col]

# Remove genes with 0 rpkms
data_c <- data_c[data_c[,tissue_name]!=0,]

# Connects to biomart to enrich our data
#listMarts()
ensembl = useMart("ensembl")
#listDatasets(ensembl)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)
# To explore the attributes and filters
# attributes[grep("gene_id", attributes[,1]),]
attributes_ensembl <- c("gene_biotype", "ensembl_gene_id", "go_id", 
                        "go_linkage_type")

# Get attributes of go_id
matrix_type <- getBM(attributes = attributes_ensembl, mart = ensembl)
matrix_type[,2:ncol(matrix_type)] <- apply(matrix_type[,2:ncol(matrix_type)], 2, as.factor)

# Get attributes of dN, dS
dN <- "mmusculus_homolog_dn"
dS <- "mmusculus_homolog_ds"
attributes_mouse_human <- c("ensembl_gene_id", dN, dS)
matrix_dn_ds <- getBM(attributes = attributes_mouse_human,  mart = ensembl)

# Merge the two datasets from Biomart
# enrich <- merge(matrix_type, matrix_dn_ds, by="ensembl_gene_id")

# Filter the dN and dS
dn_ds <- matrix_dn_ds[,dN] < 2 & matrix_dn_ds[,dN] > 0 & matrix_dn_ds[,dS] < 2 & matrix_dn_ds[,dS] > 0
dns <- dn_ds %in% TRUE
matrix_dn_ds <- matrix_dn_ds[dns,]

# Calculates the ratio
ratio <- matrix_dn_ds[,dN]/matrix_dn_ds[,dS] # Well done
ratio_matrix <- cbind(matrix_dn_ds, ratio)

# Merge data of ratios
full <- merge(data_c, ratio_matrix, by.x="Gene", by.y= "ensembl_gene_id")
full <- merge(full, matrix_type[,c("gene_biotype", "ensembl_gene_id")], by.x="Gene", by.y= "ensembl_gene_id")

#Filter by coding
proteins <- full[full$gene_biotype == "protein_coding",]
not_proteins <- full[full$gene_biotype != "protein_coding",]

res <- wilcox.test(proteins$ratio, not_proteins$ratio)
if (res$p.value < 0.01){
  print("The test is significant")
}

plot(full$gene_biotype, full$ratio, main = "The representation of ratio on type of gene",
     xlab= "Type of gene", ylab="Value of the ratio dN/dS")

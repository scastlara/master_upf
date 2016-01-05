#Script to analyse Data from project of PGB

##Change these variables accordingly
cwd <- "/home/lluis/Documents/master/PGB/exercises/Project/"
setwd(cwd)
tissue <-  "Cerebellum" # Is used to search for the tissue
# if an error occurs explore the names(raw_data) and
# store the tissue_n the number of column of your tissue
raw_file <- "GTEx_Analysis_RNA-seq_selected_tissues.txt"

#To plot some nice graphs
require(ggplot2)
# To connect to Biomart to retrieve data we can use these library 
library(biomaRt, verbose = FALSE)

# Read the data and finds the tissue column
col_classes <- c("character", "character", rep("numeric", 19))
raw_data <- read.delim(raw_file, header=TRUE, nrows = 52579, 
                     skip = 1, colClasses = col_classes)
tissue_n <-  grep(tissue, names(raw_data), ignore.case=TRUE)
tissue_name <- names(raw_data)[tissue_n]

type_gene <- function(row, tissue, specificity){
# specificity is the max number of tissue where it must be 1 or greater
# Return the tag/function of the gene
  val_all <-  sum(row >= 1)
  if (val_all == 19){
    type_g <- "Housekeeping"
  } else if (val_all <= specificity & row[tissue] >= 1){
    type_g <- "Tissue-specific"
  } else {
    type_g <- "Not-special"
  }
  return(as.vector(type_g))
}

get_genes <- function(gene){
  # Get the Ensembl identifier for each gene by spliting the version
  g <- strsplit(as.character(gene), split =".", fixed=TRUE)
  return(g[[1]][1])
}

# Select which genes are above threshold
keep <- apply(raw_data[,-c(1,2)], 1, max)
all_genes <- subset(raw_data, keep > 1)

# Remove the version of the gene
all_genes$Gene <- sapply(all_genes$Gene, get_genes)

type <- as.factor(apply(all_genes[, -c(1,2)], 1, type_gene, 
                        "Brain...Cerebellum", 3))
all_genes_type <- cbind(all_genes, type)
# Back up the data
write.table(all_genes_type, "all_genes_2.tsv", sep="\t", row.names=FALSE)


####################################

# Select the tissue genes
col <- c("Gene", "Protein", tissue_name, "type")
data_tissue <- all_genes_type[, col]

# Remove genes with 0 rpkms
data_tissue <- data_tissue[data_tissue[,tissue_name] != 0, ]
#Back up the data
write.table(data_tissue, "data_tissue_2.tsv", row.names = F)

# Connects to biomart to enrich our data
#listMarts("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
#ensembl <-  useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
#listDatasets(ensembl)
ensembl <-  useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
                    host="www.ensembl.org")
# filters <-  listFilters(ensembl)
attributes <-  listAttributes(ensembl)
# To explore the attributes and filters
#attributes[grep("gene_id", attributes[,1]),]
attributes_type <- c("ensembl_gene_id", "gene_biotype")

# Get attributes of go_id
raw_matrix_type <- getBM(attributes = attributes_type, mart = ensembl,
                      filter = "ensembl_gene_id", values = all_genes$Gene)
write.csv(raw_matrix_type, "raw_type.csv", quote=FALSE, row.names=FALSE)
# raw_matrix_type <- read.csv("raw_type.csv", header = TRUE)
matrix_type <- raw_matrix_type

protein_coding <- c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_M_gene", 
                    "IG_V_gene", "IG_Z_gene", "nonsense_mediated_decay", "nontranslating_CDS", 
                    "non_stop_decay", "polymorphic_pseudogene", "protein_coding", "TR_C_gene", 
                    "TR_D_gene", "TR_gene", "TR_J_gene", "TR_V_gene")
pseudogenes <- c("disrupted_domain", "IG_C_pseudogene", "IG_J_pseudogene", "IG_pseudogene", 
                 "IG_V_pseudogene", "processed_pseudogene", "pseudogene", 
                 "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", 
                 "translated_processed_pseudogene", "translated_unprocessed_pseudogene", 
                 "TR_J_pseudogene", "TR_V_pseudogene", "unitary_pseudogene",  "transcribed_unitary_pseudogene",
                 "unprocessed_pseudogene")
lncRNA <- c("3prime_overlapping_ncrna", "ambiguous_orf", "antisense", "lincRNA", "ncrna_host", 
            "non_coding", "processed_transcript", "retained_intron", "sense_intronic", 
            "sense_overlapping")
sncRNA <- c("miRNA", "miRNA_pseudogene", "misc_RNA", "misc_RNA_pseudogene", "Mt_rRNA", "Mt_tRNA", 
            "Mt_tRNA_pseudogene", "ncRNA", "pre_miRNA", "RNase_MRP_RNA", "RNase_P_RNA", "rRNA", 
            "rRNA_pseudogene", "scRNA_pseudogene", "snlRNA", "snoRNA", "snoRNA_pseudogene", "snRNA", 
            "snRNA_pseudogene", "SRP_RNA", "tmRNA", "tRNA", "tRNA_pseudogene", "scaRNA")

# Add the type of gene to the data
full_tissue <- merge(data_tissue, matrix_type, by.x="Gene", by.y= "ensembl_gene_id",
              all.y = FALSE, all.x=TRUE, sort = FALSE)

full_all <- merge(all_genes_type, matrix_type, by.x="Gene", by.y= "ensembl_gene_id",
                     all.y = FALSE, all.x=TRUE, sort = FALSE)
# Reclasify the biotypes
full_tissue$gene_biotype <- as.character(full_tissue$gene_biotype)
full_tissue[full_tissue$gene_biotype %in% protein_coding, "gene_biotype"] <- "protein_coding"
full_tissue[full_tissue$gene_biotype %in% lncRNA, "gene_biotype"] <- "lncRNA"
full_tissue[full_tissue$gene_biotype %in% pseudogenes, "gene_biotype"] <- "pseudogenes"
full_tissue[full_tissue$gene_biotype %in% sncRNA, "gene_biotype"] <- "sncRNA"
full_tissue[full_tissue$gene_biotype %in% c("TEC", NA), "gene_biotype"] <- "unknown"
full_tissue$gene_biotype <- as.factor(full_tissue$gene_biotype)

full_all$gene_biotype <- as.character(full_all$gene_biotype)
full_all[full_all$gene_biotype %in% protein_coding, "gene_biotype"] <- "protein_coding"
full_all[full_all$gene_biotype %in% lncRNA, "gene_biotype"] <- "lncRNA"
full_all[full_all$gene_biotype %in% pseudogenes, "gene_biotype"] <- "pseudogenes"
full_all[full_all$gene_biotype %in% sncRNA, "gene_biotype"] <- "sncRNA"
full_all[full_all$gene_biotype %in% c("TEC", NA), "gene_biotype"] <- "unknown"
full_all$gene_biotype <- as.factor(full_all$gene_biotype)

# Get attributes of dN, dS
attributes_mouse_human <- c("ensembl_gene_id", "mmusculus_homolog_dn",
                            "mmusculus_homolog_ds")
raw_matrix_dn_ds <- getBM(attributes = attributes_mouse_human,  mart = ensembl,
                      filter = "ensembl_gene_id", values = all_genes$Gene)
write.csv(raw_matrix_dn_ds, "raw_dn_ds.csv", quote=FALSE, row.names=FALSE)
# raw_matrix_dn_ds <- read.csv("raw_dn_ds.csv")
# Filter the dN and dS
dn_ds <- raw_matrix_dn_ds$mmusculus_homolog_dn < 2 & raw_matrix_dn_ds$mmusculus_homolog_ds < 2
dns <- dn_ds %in% TRUE
matrix_d <- raw_matrix_dn_ds[dns,]

# Calculates the ratio
ratio <- matrix_d$mmusculus_homolog_dn/matrix_d$mmusculus_homolog_ds # Well done
ratio_matrix <- cbind(matrix_d, ratio)
keep <- ratio_matrix$ratio < 2
keep <- keep %in% TRUE
ratio_matrix <- ratio_matrix[keep,]

# If a same gene Id has two ratio order them
ordered_marix_ratio <- ratio_matrix[order(ratio_matrix$ensembl_gene_id,
                                          -rank(ratio_matrix$ratio), decreasing = TRUE),]
# Select the lower score of the ratio
right_dn_ds <- ordered_marix_ratio[
  match(unique(ordered_marix_ratio$ensembl_gene_id),
        ordered_marix_ratio$ensembl_gene_id),]

# Add the dn_ds to the type
full_tissue <- merge(full_tissue, right_dn_ds, by.x="Gene", by.y= "ensembl_gene_id",
              all.y = FALSE, all.x=TRUE, sort = FALSE)
full_all <- merge(full_all, right_dn_ds, by.x="Gene", by.y= "ensembl_gene_id",
                     all.y = FALSE, all.x=TRUE, sort = FALSE)

# Back up the data without GO but with dS and dS
write.table(full_tissue, "full_tissue_2.tsv", sep="\t", row.names=FALSE)
write.table(full_all, "full_all_2.tsv", sep="\t", row.names=FALSE)

# Save the Data loaded
save.image("Project_2.RData")

################################################
## Evolution questions/relationships

pdf("images_tissue_2.pdf")
# postscript()
# General configuration for ggplot
p <- ggplot(full_tissue) # Genes in our tissue expressing
pt <- ggplot(full_all) # All genes
pt_t <- ggplot(full_all[full_all$type == "Tissue-specific",]) #All genes tissue-specific
p_t <- ggplot(full_tissue[full_tissue$type == "Tissue-specific",]) # Expressing genes tissue-specific
l <- theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 1)) # Options of the labels
b <-  theme(axis.line = element_line(colour = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(colour = "grey", size = 0.5, linetype = 2),
            panel.grid.minor.y = element_line(colour = "grey", size = 0.25, linetype = 3),
            panel.border = element_blank(),
            panel.background = element_blank()) 

# Check if there is non normality
p + l + b +
  geom_bar(aes(gene_biotype, fill=type)) + 
  ggtitle("Accumulative biotype for each type of gene in our tissue") +
  labs(x = "Biotype")
pt + l + b +
  geom_bar(aes(gene_biotype, fill=type)) + 
  ggtitle("Accumulative biotype for all genes") +
  labs(x = "Biotype") 
p + b +
  geom_bar(aes(gene_biotype, fill=type), position="fill") +
  ggtitle("Fraction of each biotype on each type of gene with expression in our tissue") +
  labs(x="Biotype", y="Frequency")

pt + b +
  geom_bar(aes(gene_biotype, fill=type), position="fill") +
  ggtitle("Fraction of each biotype on each type of gene in our tissue") +
  labs(x="Biotype", y="Frequency") 
p + b +
  geom_bar(aes(type, fill=gene_biotype), position="fill") +
  ggtitle("Fraction of each biotype on each type of expressing gene in our tissue") +
  labs(x="Type of gene", y="Frequency")
pt + b +
  geom_bar(aes(type, fill=gene_biotype), position="fill") +
  ggtitle("Fraction of each biotype on each type of gene") +
  labs(x="Type of gene", y="Frequency")
p + l + b + 
  geom_bar(aes(gene_biotype)) + 
  facet_wrap(~ type) + 
  ggtitle("Number biotype for each type of expressing gene in our tissue") +
  labs(x = "Biotype")
pt + l + b +
  geom_bar(aes(gene_biotype)) + 
  facet_wrap(~ type) + 
  ggtitle("Number biotype for each type of gene") +
  labs(x = "Biotype")
p + l + b +
  geom_bar(aes(type)) + 
  facet_grid(~ gene_biotype) + 
  ggtitle("Number of type of expressing gene for each biotype") +
  labs(x = "Type of gene")
pt + l + b +
  geom_bar(aes(type)) + 
  facet_grid(~ gene_biotype) + 
  ggtitle("Number of type of gene for each biotype") +
  labs(x = "Type of gene") 
qplot(sample = full_tissue$mmusculus_homolog_dn) + 
  ggtitle("Distribution of dN") +
  b
qplot(sample = full_tissue$mmusculus_homolog_ds) + 
  ggtitle("Distribution of dS") +
  b
qqplot(full_tissue$mmusculus_homolog_dn, full_tissue$mmusculus_homolog_ds,
       main="QQ plot of the representation of dS over dN",
       xlab = "Value of dN", ylab = "Value of dS", col = full_tissue$type)
p + b +
  geom_boxplot(aes(type,ratio)) + 
  ggtitle("Evolutionary pressure over type of genes") + 
  labs(x="Type of gene", y = "Ratio dN/dS") +
  ylim(0, 2) 
# p + b + l +
#   geom_boxplot(aes(type,ratio)) + 
#   ggtitle("Evolutionary pressure over type of genes") + 
#   labs(x="Type of gene", y = "Ratio dN/dS") +
#   facet_wrap(~ gene_biotype)+
#   ylim(0, 2)
# pt + b + l +
#   geom_boxplot(aes(type,ratio)) + 
#   ggtitle("Evolutionary pressure over type of genes") + 
#   labs(x="Type of gene", y = "Ratio dN/dS") +
#   facet_wrap(~ gene_biotype)+
#   ylim(0, 2)
# plot(full_tissue$type, full_tissue$ratio,
#      main= "Does the evulutionary pressure affect differently the genes",
#      xlab="Type of genes", ylab="Value of the dN/dS")

# Ploting data bout the evolutionary question
# summary(full_tissue$gene_biotype[!is.na(full_tissue$ratio)])
p + b +
  geom_boxplot(aes(gene_biotype, ratio)) +
  ggtitle("Ratio for type of gene") +
  labs(x= "Biotype", y = "Ratio dN/dS") +
  ylim(0, 2)
# op <- par(mar=c(14,4,4,5)+0.1)
# plot(full_tissue$gene_biotype, full_tissue$ratio, main = "The representation of ratio on type of gene",
#      ylab = "Value of the ratio dN/dS")
# par(mar=c(5,4,4,2)+0.1)
pt + b +
  geom_point(aes(mmusculus_homolog_ds, mmusculus_homolog_dn, colour=factor(type))) + 
  ggtitle("dN over dS") +
  labs(x= "dS", y = "dN") +
  geom_abline(intercept = 0, slope=1)
plot(full_tissue$mmusculus_homolog_ds, full_tissue$mmusculus_homolog_dn, main="Representation of dN/dS",
      xlab = "Value of dS", ylab = "Value of dN", col = full_tissue$type) 
abline(0,1)
# legend("topleft", legend = c("House keeping", "Tissue-specific", "Not relevant"),
#        col=1:length(data$type),pch=1, cex=0.6, xjust = 1, yjust = 1)

# The distribution of each sample
hist(full_tissue$mmusculus_homolog_dn, main = "Histogram of dN", xlab="dN", ylab="Counts")
hist(full_tissue$mmusculus_homolog_ds, main = "Histogram of dS", xlab="dS", ylab="Counts")
# op <- par(mar=c(15,4,4,2))
# pdf("test.pdf")
# barplot(prop.table(table(full_tissue[(is.na(full_tissue$ratio)),"gene_biotype"])),
#         main = "Biotype distributed without ratio dN/dS",
#         ylab = "Percentatge of genes",
#         las = 2)
# dev.off()
# op


########
# Calculus
# Linear correlation dN-dS
cor(full_tissue$mmusculus_homolog_dn[!is.na(full_tissue$mmusculus_homolog_dn)],
    full_tissue$mmusculus_homolog_ds[!is.na(full_tissue$mmusculus_homolog_ds)])
kruskal.test(mmusculus_homolog_dn ~ mmusculus_homolog_ds, data = full_tissue) 
#Significative
# Anova test of which factors affect the dN/dS
anova(lm(ratio ~ type, data = full_tissue))
#significative
kruskal.test(ratio ~ type, data = full_tissue[full_tissue$type == "Housekeeping",]) 
# Same group
kruskal.test(ratio ~ gene_biotype, data = full_tissue[full_tissue$type == "Housekeeping",]) 
#Same group
kruskal.test(ratio ~ Brain...Cerebellum, data = full_tissue[full_tissue$type == "Housekeeping",]) 
# p-value 0.4854
anova(lm(ratio ~ gene_biotype, data = full_tissue))
anova(lm(ratio ~ gene_biotype*type, data = full_tissue))

# As we can see it is not a straight line then we must use the wilcox test
#Filter by coding
proteins <- full_tissue[full_tissue$gene_biotype == "protein_coding",]
not_proteins <- full_tissue[full_tissue$gene_biotype != "protein_coding",]
# Compare if the ratio is greater in tissue-specific genes against others
wilcox.test(proteins[full_tissue$type == "Tissue-specific", "ratio"],
            proteins[full_tissue$type != "Tissue-specific", "ratio"], alternative = "greater")
# p-value 0.9453
# DEG
# Transcript més llarg biomart ??
# Normalitzat per isoformas i summat. 
anova(lm(full_tissue$Brain...Cerebellum~full_tissue$gene_biotype))
#1.01e-12
kruskal.test(Brain...Cerebellum~gene_biotype,data=full_tissue)
#p-value 2.2e-16
anova(lm(full_tissue$Brain...Cerebellum~full_tissue$type))
#2.2e-16
kruskal.test(Brain...Cerebellum~type,data=full_tissue)
#2.2e-16
anova(lm(full_tissue$Brain...Cerebellum~full_tissue$type*full_tissue$gene_biotype))
# All significant
foldc <- full_tissue[,"Brain...Cerebellum"]/mean(
  full_tissue[full_tissue$type == "Housekeeping", "Brain...Cerebellum"])
full_foldc<- cbind(full_tissue, foldc)
write.csv(full_foldc, "full_tissue_fc.csv", row.names=FALSE, quote=FALSE)
f <- ggplot(full_foldc)
f + b +
  geom_point(aes(type, log10(foldc), color = gene_biotype))
f + b +
  geom_boxplot(aes(type, log10(foldc))) +
  ggtitle("Comparing the mean expression of the housekeeping with the genes") +
  labs(x="Gene type")
# Proportion of data explained by the type, the gene_biotype, the interaction and the residual variation
prop.table(multi$`Sum Sq`) 
# compare proteins and non proteins
wilcox.test(proteins$Brain...Cerebellum,
            not_proteins$Brain...Cerebellum, "greater")
#2.2e-16
t_proteins <- proteins[proteins$type == "Tissue-specific", ]
not_t_proteins <- proteins[proteins$type != "Tissue-specific", ]
wilcox.test(t_proteins$Brain...Cerebellum,
            not_t_proteins$Brain...Cerebellum, "less")
#0.0003957
wilcox.test(full_all[full_all$gene_biotype=="sncRNA" & full_all$type == "Tissue-specific", "Brain...Cerebellum"],
            full_all[full_all$gene_biotype!="sncRNA"& full_all$type == "Tissue-specific", "Brain...Cerebellum"], "less")
# 0.0001376
# Ploting respect the level of expression of the data
p + 
  geom_point(aes(log10(Brain...Cerebellum), ratio, color = type)) + 
  ggtitle("Evolution pressure over expression") +
  xlab("log10(rpkm)")+
  b
p + 
  geom_point(aes(log10(Brain...Cerebellum), ratio, color = gene_biotype)) + 
  ggtitle("Evolution pressure over expression") +
  xlab("log10(rpkm)")+
  l +
  b

p_t + 
  geom_point(aes(log10(Brain...Cerebellum), ratio, color = type)) + 
  ggtitle("Tissue-specific genes") +
  xlab("log10(rpkm)")+
  l + b

# plot(log10(full_all$Brain...Cerebellum), full_all$ratio)
#TODO see if this representation explain if protein_coding are more or less expressed than others
p +
  geom_boxplot(aes(gene_biotype, log10(Brain...Cerebellum))) +
  ggtitle("Expression in cerebellum by biotype") +
  labs(x="Gene Biotype", y= "log10(rpkm)") +
  l+
  b

p_t +
  geom_boxplot(aes(gene_biotype, log10(Brain...Cerebellum))) +
  ggtitle("Expression in cerebellum by biotype in tissue-specific genes") +
  labs(x="Gene Biotype", y= "log10(rpkm)") +
  b

qqnorm(log10(proteins$Brain...Cerebellum), c(-4,4),
       main = "Normal Q-Q plot of the log10 of rpkm of proteins")
qqline(log10(proteins$Brain...Cerebellum))

qqnorm(log10(not_proteins$Brain...Cerebellum), c(-4,4),
       main = "Normal Q-Q plot of the log10 of rpkm of not proteins")
qqline(log10(not_proteins$Brain...Cerebellum))

p + 
  geom_boxplot(aes(gene_biotype, log10(Brain...Cerebellum))) + 
  facet_grid(~type) +
  ggtitle("Expression in Cerebellum") + 
  labs(x="Biotype", y = "log10(rpkm)") +
  l + 
  b

p + 
  geom_boxplot(aes(type, log10(Brain...Cerebellum))) + 
  facet_grid(~gene_biotype) +
  ggtitle("Expression in Cerebellum") + 
  labs(x="Types of gene", y = "log10(rpkm)") +
  l + 
  b

# plot(full_tissue$gene_biotype, log10(full_tissue$Brain...Cerebellum), ylab="log10(rpkm)",
#      las=2, main = "log10arithm of rpkm respect gene biotype in Cerebellum")

# plot(full_tissue$type, log10(full_tissue$Brain...Cerebellum), 
#      main = "log10arithmic values of rpkm for each type of tissue",
#      xlab = "Tissue", ylab = "log10(rpkm)")
freq <- sort(prop.table(table(full_tissue[full_tissue$type == "Tissue-specific", "gene_biotype"])),
             decreasing = T)
pie(freq, main = "Tissue specific representation")
freq <- sort(prop.table(table(full_tissue[full_tissue$type == "Housekeeping", "gene_biotype"])),
             decreasing = T)
pie(freq, main = "Housekeeping representation")
freq <- sort(prop.table(table(full_tissue[full_tissue$type == "Not-special", "gene_biotype"])),
             decreasing = T)
pie(freq, radius=1, main = "Not-special representation")
freq <- sort(prop.table(table(full_tissue[full_tissue$type != "Tissue-specific", "gene_biotype"])))
pie(freq, main = "Not tissue-specific representation")
# legend("topleft", legend = c("Protein coding", "lncRNA", "snRNA", "pseudogenes"))
dev.off()

# GO
# Analysis of the Gene Ontology 

#Get attributes of GO
attributes_go <- c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003")
#"definition_1006", "go_linkage_type", "namespace_1003",
#"goslim_goa_accession", "goslim_goa_description",
#"go_molecular_function_id", "go_molecular_function__dm_name_1006")

# Load data from biomart
raw_matrix_go <- getBM(attributes = attributes_go, mart = ensembl,
                       filter = "ensembl_gene_id", values = all_genes$Gene)
write.csv(raw_matrix_go, "raw_go.csv", quote=FALSE, row.names=FALSE)
# raw_matrix_go <- read.csv("raw_go.csv")
# Remove empty data
matrix_go <- raw_matrix_go[raw_matrix_go$go_id != "",]
# Select tissue-specific
ts_genes <- full_tissue[full_tissue$type == "Tissue-specific", "Gene"]
ts_p_genes <- full_tissue[full_tissue$type == "Tissue-specific" & full_tissue$gene_biotype == "protein_coding", "Gene"]

# Calculate the go_id frequency  of the intersting genes
matrix_go_1 <- subset(matrix_go, unique(matrix_go$ensembl_gene_id) %in% ts_genes)
matrix_go_1$namespace_1003 <- as.factor(matrix_go_1$namespace_1003)


matrix_go_2 <- subset(matrix_go, unique(matrix_go$ensembl_gene_id) %in% ts_p_genes)
matrix_go_2$namespace_1003 <- as.factor(matrix_go_2$namespace_1003)

# Define function to calculate frequency of each GO
select_bioprocess <- function(matrix, process){
  # select those specific type
  t <- matrix[matrix$namespace_1003 == process, ] 
  # Calculate the frequency and the acumulative frequency
  freq_t <- sort(prop.table(table(t$go_id)), T)
  ac_freq_t <- cumsum(freq_t)
  #Create a data.frame with these data
  go_freq_t <- data.frame(go_id =names(freq_t), 
                           frequency = freq_t, 
                           acumulative_freq = ac_freq_t)
  # merge the data to link frequency and descriptions
  go_descr <- merge(go_freq_t, t[,-c(1,4)], by.x = "go_id", 
                    by.y = "go_id", all.x = TRUE, all.y = FALSE, sort = FALSE)
  go_done <- go_descr[!duplicated(go_descr),]
  go_done$frequency <- as.numeric(go_done$frequency)
  go_done$acumulative_freq <- as.numeric(go_done$acumulative_freq)
  return(go_done)
}

bp <- select_bioprocess(matrix_go_1, "biological_process")
cc <- select_bioprocess(matrix_go_1, "cellular_component")
mf <- select_bioprocess(matrix_go_1, "molecular_function")
write.csv(bp, "go_bp.csv", row.names = FALSE)
write.csv(cc, "go_cc.csv", row.names = FALSE)
write.csv(mf, "go_mf.csv", row.names = FALSE)

bp <- select_bioprocess(matrix_go_2, "biological_process")
cc <- select_bioprocess(matrix_go_2, "cellular_component")
mf <- select_bioprocess(matrix_go_2, "molecular_function")
write.csv(bp, "go_protein_bp.csv", row.names = FALSE)
write.csv(cc, "go_protein_cc.csv", row.names = FALSE)
write.csv(mf, "go_protein_mf.csv", row.names = FALSE)
save.image("Project_2.RData")

############################
# Length affects it

raw_matrix_length <- getBM(attributes = c("ensembl_gene_id","transcript_length"), mart = ensembl,
                           filter = "ensembl_gene_id", values = all_genes$Gene)
write.csv(raw_matrix_length, "raw_length.csv", quote=FALSE, row.names=FALSE)
# raw_matrix_length <- read.csv("raw_length.csv")
matrix_length <- aggregate( raw_matrix_length$transcript_length ~ raw_matrix_length$ensembl_gene_id,
           raw_matrix_length, mean)
names(matrix_length) <- c("ensembl_gene_id", "length")
full_all_length <- merge(full_all, matrix_length, 
                         by.x= "Gene", by.y = "ensembl_gene_id",
                         all.x = TRUE, all.y=FALSE)
full_all_length$gene_biotype <- as.factor(full_all_length$gene_biotype)
attach(full_all_length)
# ANCOVA
analysis <- anova(lm(full_all_length$Brain...Cerebellum ~ full_all_length$length+full_all_length$type*full_all_length$gene_biotype))
analysis_2 <- anova(lm(full_all_length$Brain...Cerebellum ~ full_all_length$length*full_all_length$type*full_all_length$gene_biotype))
wilcox.test(full_all_length$Brain...Cerebellum, full_all_length$length, paired=TRUE, "less")
#2.2e-16
wilcox.test(full_all_length$ratio, full_all_length$length, paired=TRUE)
#2.2e-16
kruskal.test(full_all_length$Brain...Cerebellum~full_all_length$type)
#2.2e-16
kruskal.test(full_all_length$length~full_all_length$type, paired=TRUE, "less")
kruskal.test(full_all_length$length~full_all_length$gene_biotype, paired=TRUE, "less")
kruskal.test(full_all_length$length~full_all_length$type, paired=TRUE)
anova(lm(full_all_length$length ~ full_all_length$type*full_all_length$gene_biotype))

len <- ggplot(full_all_length)
pdf("rpkm_length.pdf")
len +
  geom_point(aes(log10(length), log10(Brain...Cerebellum), col=type)) +
  ggtitle("Expression of genes over the length") +
  ylab("Log10(rpkm)") +
  b

len +
  geom_point(aes(log10(length), log10(Brain...Cerebellum), col=gene_biotype)) +
  ggtitle("Expression of genes over the length") +
  ylab("Log10(rpkm)") +
  b

len_t <- ggplot(full_all_length[full_all_length$type == "Tissue-specific",])
len_t + 
  geom_point(aes(log10(length), log10(Brain...Cerebellum), col=gene_biotype)) +
  ggtitle("Tissue-specific genes") +
  xlab("log10(length)")+
  ylab("log10(rpkm)")+
  b
dev.off()

save.image("Project_2.RData")

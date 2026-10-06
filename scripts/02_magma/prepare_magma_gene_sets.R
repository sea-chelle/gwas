# prepare experimental gene sets for MAGMA gene-set analysis
#
# gene sets are derived from RNA-seq differential expression,
# WGCNA and ATAC-seq/PLAC-seq analyses in mouse hypothalamus.
# mouse genes are mapped to high-confidence human orthologues
# and converted to Entrez IDs represented in the MAGMA NCBI37.3
# gene reference.

# get required packages
library(tidyverse)
library(biomaRt)

# input files
wgcna_file <- "path/to/WGCNA_grey60.rds"
magma_gene_file <- "path/to/NCBI37.3_gene_ids.tsv"

# output file
output_file <- "stress_gene_sets.set.annot"

# create RNA-seq gene sets from significant DEGs
rna_gene_sets_mouse <- list(RNA_Overall_60 = rna_results$basal_vs_60 %>%
    filter(diffexpressed != "NO") %>%
    pull(ensembl) %>%
    unique(),
  
  RNA_Male_60 = rna_results$male_basal_vs_60 %>%
    filter(diffexpressed != "NO") %>%
    pull(ensembl) %>%
    unique(),
  
  RNA_Female_60 = rna_results$female_basal_vs_60 %>%
    filter(diffexpressed != "NO") %>%
    pull(ensembl) %>%
    unique())

# check number of genes in each RNA-seq gene set
sapply(rna_gene_sets_mouse, length)

# import WGCNA grey60 module genes
WGCNA_grey60 <- readRDS(wgcna_file)

wgcna_gene_set_mouse <- unique(WGCNA_grey60$gene)

# check number of genes in WGCNA gene set
length(wgcna_gene_set_mouse)

# create ATAC-seq gene sets from PLAC-seq-linked promoter genes
atac_gene_sets_mouse <- list(ATAC_Overall_60 = comparisons$basal_vs_60$plac_peak_gene %>%
                               pull(ENSEMBL_ID) %>%
                               na.omit() %>%
                               unique(),
  
                             ATAC_Male_60 = comparisons$male_basal_vs_60$plac_peak_gene %>%
                               pull(ENSEMBL_ID) %>%
                               na.omit() %>%
                               unique(),

                             ATAC_Female_60 = comparisons$female_basal_vs_60$plac_peak_gene %>%
                               pull(ENSEMBL_ID) %>%
                               na.omit() %>%
                               unique())

# check number of genes in each ATAC-seq gene set
sapply(atac_gene_sets_mouse, length)

# combine RNA-seq, WGCNA and ATAC-seq gene sets
gene_sets_mouse <- c(rna_gene_sets_mouse,
                     list(WGCNA_grey60 = wgcna_gene_set_mouse),
                     atac_gene_sets_mouse)

# check number of genes in each mouse gene set
sapply(gene_sets_mouse, length)

# get all unique mouse genes across gene sets
all_mouse_genes <- unique(unlist(gene_sets_mouse))

# connect to Ensembl mouse dataset
mouse_mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# map mouse Ensembl genes to human orthologues
orthologue_map <- getBM(attributes = c("ensembl_gene_id",
                                       "external_gene_name",
                                       "hsapiens_homolog_ensembl_gene",
                                       "hsapiens_homolog_associated_gene_name",
                                       "hsapiens_homolog_orthology_type",
                                       "hsapiens_homolog_orthology_confidence"),
                        filters = "ensembl_gene_id",
                        values = all_mouse_genes,
                        mart = mouse_mart)

# simplify column names
colnames(orthologue_map) <- c("mouse_ensembl", 
                              "mouse_symbol", 
                              "human_ensembl", 
                              "human_symbol", 
                              "orthology_type", 
                              "orthology_confidence")

# retain high-confidence human orthologues
orthologue_map_highconf <- orthologue_map %>%
  filter(!is.na(orthology_confidence),
         orthology_confidence == 1,
         human_ensembl != "") %>%
  distinct()

# check orthology types and mapping success
table(orthologue_map_highconf$orthology_type)
length(unique(orthologue_map_highconf$mouse_ensembl))
length(unique(orthologue_map_highconf$human_ensembl))

# calculate percentage of mouse genes with a high-confidence human orthologue
n_mouse_input <- length(all_mouse_genes)
n_mouse_mapped <- length(unique(orthologue_map_highconf$mouse_ensembl))
100 * n_mouse_mapped / n_mouse_input

# convert mouse gene sets to human Ensembl gene sets
gene_sets_human <- lapply(gene_sets_mouse, function(mouse_genes) {
  orthologue_map_highconf %>%
    filter(mouse_ensembl %in% mouse_genes) %>%
    pull(human_ensembl) %>%
    unique()
})

# compare mouse and human gene set sizes
mapping_summary <- data.frame(gene_set = names(gene_sets_mouse),
                              mouse_n = sapply(gene_sets_mouse, length),
                              human_n = sapply(gene_sets_human, length))

mapping_summary

# check human Ensembl gene IDs
lapply(gene_sets_human, head)

# get all unique human genes across gene sets
all_human_genes <- unique(unlist(gene_sets_human))

# connect to Ensembl human dataset
human_mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# convert human Ensembl IDs to Entrez IDs
human_id_map <- getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name",
                                     "entrezgene_id"),
                      filters = "ensembl_gene_id", values = all_human_genes, mart = human_mart)

# simplify column names
colnames(human_id_map) <- c("human_ensembl", "human_symbol", "entrez_id")

# import genes represented in the MAGMA NCBI37.3 reference
magma_genes <- read.delim(magma_gene_file, header = FALSE, col.names = c("entrez_id", "magma_symbol"))

# convert Entrez IDs to character
human_id_map <- human_id_map %>%
  mutate(entrez_id = as.character(entrez_id))

magma_genes <- magma_genes %>%
  mutate(entrez_id = as.character(entrez_id))

# retain human genes represented in the MAGMA NCBI37.3 reference
human_id_map_magma <- human_id_map %>%
  filter(!is.na(entrez_id),
         entrez_id %in% magma_genes$entrez_id) %>%
  distinct()

# check number of human genes retained in the MAGMA reference
length(unique(human_id_map$human_ensembl))
length(unique(human_id_map_magma$human_ensembl))

length(unique(na.omit(human_id_map$entrez_id)))
length(unique(human_id_map_magma$entrez_id))

# convert human gene sets to MAGMA-compatible Entrez IDs
gene_sets_magma <- lapply(gene_sets_human, function(human_genes) {
  human_id_map_magma %>%
    filter(human_ensembl %in% human_genes) %>%
    pull(entrez_id) %>%
    unique()
})

# check final MAGMA gene set sizes
sapply(gene_sets_magma, length)

# check MAGMA gene sets contain Entrez IDs
lapply(gene_sets_magma, head)

# write gene sets in MAGMA set-annotation format
writeLines(vapply(names(gene_sets_magma), function(set_name) {
      paste(set_name,
            paste(gene_sets_magma[[set_name]], collapse = "\t"),
            sep = "\t")
    },
    character(1)),
    output_file)

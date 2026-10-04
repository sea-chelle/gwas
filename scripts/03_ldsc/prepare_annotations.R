# prepare experimental annotations for S-LDSC analysis
#
# RNA-seq annotations are derived from significant DEGs in mouse hypothalamus.
# mouse genes are mapped to high-confidence human orthologues and converted to
# GRCh37 genomic regions spanning the gene body +/- 100 kb.
#
# ATAC-seq annotations are derived directly from significant DAR coordinates
# and converted to BED format for cross-species liftOver from mm39 to hg19.

# get required packages
library(tidyverse)
library(biomaRt)

# set working directory
setwd(proj_dir)

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

# get all unique mouse genes across RNA-seq gene sets
all_mouse_genes <- unique(unlist(rna_gene_sets_mouse))

# connect to Ensembl mouse dataset
mouse_mart <- useEnsembl(biomart = "genes",
                         dataset = "mmusculus_gene_ensembl")

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

# convert mouse RNA-seq gene sets to human Ensembl gene sets
gene_sets_human <- lapply(rna_gene_sets_mouse, function(mouse_genes) {
  orthologue_map_highconf %>%
    filter(mouse_ensembl %in% mouse_genes) %>%
    pull(human_ensembl) %>%
    unique()
})

# check number of genes in each human gene set
sapply(gene_sets_human, length)

# get all unique human genes across gene sets
all_human_genes <- unique(unlist(gene_sets_human))

# connect to GRCh37 Ensembl human dataset
human_mart_grch37 <- useEnsembl(biomart = "genes",
                                dataset = "hsapiens_gene_ensembl",
                                GRCh = 37)

# retrieve GRCh37 genomic coordinates
human_gene_coords <- getBM(attributes = c("ensembl_gene_id",
                                          "chromosome_name",
                                          "start_position",
                                          "end_position"),
                           filters = "ensembl_gene_id",
                           values = all_human_genes,
                           mart = human_mart_grch37)

# check output
dim(human_gene_coords)
head(human_gene_coords)

# check how many human genes there were originally
length(all_human_genes)

# check how many human genes were mapped to GRCh37 coordinates
length(unique(human_gene_coords$ensembl_gene_id))

# check how many genes are on each chromosome
table(human_gene_coords$chromosome_name)

# retain autosomal genes represented by the LDSC reference
human_gene_coords_autosomal <- human_gene_coords %>%
  filter(chromosome_name %in% as.character(1:22)) %>%
  distinct()

# check retained genes
nrow(human_gene_coords_autosomal)
length(unique(human_gene_coords_autosomal$ensembl_gene_id))

# identify genes excluded because they are non-autosomal
human_gene_coords %>%
  filter(!chromosome_name %in% as.character(1:22))

# check how many genes in each set remain after mapping and filtering non-autosomal genes
sapply(gene_sets_human, function(genes) {
  sum(genes %in% human_gene_coords_autosomal$ensembl_gene_id)
})

# create +/- 100 kb GRCh37 regions for each RNA-seq gene set
rna_gene_regions <- lapply(gene_sets_human, function(genes) {
  
  human_gene_coords_autosomal %>%
    filter(ensembl_gene_id %in% genes) %>%
    mutate(start = pmax(1, start_position - 100000),
           end = end_position + 100000) %>%
    dplyr::select(ensembl_gene_id,
                  chromosome = chromosome_name,
                  start,
                  end) %>%
    arrange(as.numeric(chromosome), start)
})

# convert RNA-seq regions from 1-based genomic coordinates to 0-based BED format
rna_bed <- lapply(rna_gene_regions, function(x) {
  x %>%
    transmute(chr = paste0("chr", chromosome),
              start = start - 1L,
              end = end,
              gene = ensembl_gene_id)
})

# check BED-formatted RNA-seq regions
sapply(rna_bed, nrow)
head(rna_bed$RNA_Overall_60)

# prepare ATAC-seq DAR coordinates for liftOver
atac_regions_mouse <- lapply(master_dar_tables, function(x) {
  x %>%
    dplyr::select(seqnames, start, end, peak) %>%
    distinct()
})

# check number of DARs in each ATAC-seq annotation
sapply(atac_regions_mouse, nrow)
lapply(atac_regions_mouse, head)

# convert DARs from 1-based genomic coordinates to 0-based BED format for liftOver
atac_bed <- lapply(master_dar_tables, function(x) {
  x %>%
    transmute(chr = as.character(seqnames),
              start = start - 1L,
              end = end,
              peak = peak) %>%
    distinct()
})

# check BED-formatted DARs
sapply(atac_bed, nrow)
head(atac_bed$Overall)

# create output directory
annotation_dir <- file.path(proj_dir, "ldsc/annotations_input")
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)

# export RNA-seq GRCh37 regions as BED4 files
for (set_name in names(rna_bed)) {
  
  write.table(rna_bed[[set_name]],
              file = file.path(annotation_dir,
                               paste0(set_name, "_GRCh37_100kb.bed")),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}

# export ATAC-seq mm39 regions as BED4 files
for (group in names(atac_bed)) {
  
  write.table(atac_bed[[group]],
              file = file.path(annotation_dir,
                               paste0("ATAC_", group, "_mm39.bed")),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}

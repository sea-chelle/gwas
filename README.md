# GWAS integration with transcriptomic and epigenomic data

This repository contains analysis workflows for integrating psychiatric genome-wide association study (GWAS) summary statistics with transcriptomic and epigenomic datasets using MAGMA and stratified LD score regression (S-LDSC). The pipeline supports gene-set and genomic annotation-based analyses for investigating relationships between GWAS genetic associations and features derived from RNA-seq and ATAC-seq data.

## Analyses

### MAGMA
MAGMA is used to perform gene-level association analysis and test enrichment of GWAS association within user-defined gene sets. Gene sets can be derived from transcriptomic data or regulatory regions linked to putative target genes.

### S-LDSC
S-LDSC is used to test whether genomic annotations derived from transcriptomic or epigenomic analyses explain disproportionate SNP heritability. Custom annotations are analysed conditional on the baselineLD v2.3 model. The workflow also supports joint S-LDSC models for directly comparing related annotations using regression coefficient covariance. Examples include sex (male vs female comparisons).

## Data and reference resources
Raw GWAS summary statistics and large reference datasets are not distributed with this repository and must be downloaded separately.

### GWAS summary statistics
The current manifest includes summary statistics for:
- Attention-deficit/hyperactivity disorder (ADHD)
- Anxiety disorders (ANX)
- Autism spectrum disorder (ASD)
- Bipolar disorder (BD)
- Major depressive disorder (MDD)
- Post-traumatic stress disorder (PTSD)
- Schizophrenia (SCZ)

Details of the GWAS study, expected input file, genome build, column names and sample-size handling are provided in `manifests/GWAS_manifest.tsv`. Raw GWAS summary statistics are not redistributed in this repository. Summary statistics can be obtained from the [Psychiatric Genomics Consortium (PGC)](https://pgc.unc.edu/for-researchers/download-results/).

### MAGMA reference data
MAGMA requires the MAGMA executable, gene-location files and an appropriate genotype reference panel for gene analysis.

This workflow uses:
- MAGMA v1.10
- NCBI37.3 gene locations
- European 1000 Genomes Phase 3 reference genotypes

These resources can be obtained from the [MAGMA website](https://cncr.nl/research/magma/).

### S-LDSC reference data
S-LDSC analyses require reference genotype and LD-score resources. This workflow uses:
- European 1000 Genomes Phase 3 PLINK files
- European 1000 Genomes Phase 3 allele frequencies
- HapMap3 SNPs excluding the MHC
- European regression weights
- baselineLD v2.3 annotations and LD scores

These reference files are available from the [LDSC reference-data repository.](https://zenodo.org/records/10515792). S-LDSC analyses in this workflow use the [CBIIT Python 3 implementation of LDSC](https://github.com/CBIIT/ldsc).

## Repository structure

```text
gwas/
├── README.md
├── environment/
│   └── environment.yml
├── manifests/
│   └── GWAS_manifest.tsv
├── scripts/
│   ├── 01_prepare_gwas/
│   │   ├── prepare_asd.sh
│   │   └── prepare_gwas.sh
│   ├── 02_magma/
│   │   ├── annotate_genes.sh
│   │   ├── gene_analysis.sh
│   │   └── gene_set_analysis.sh
│   ├── 03_ldsc/
│   │   ├── prepare_annotations.R
│   │   ├── liftover_atac.sh
│   │   ├── prepare_gwas_ldsc.sh
│   │   ├── run_sldsc.sh
│   │   ├── run_sldsc_array.sh
│   │   ├── run_sldsc_sex_comparison.sh
│   │   └── run_sldsc_sex_comparison_array.sh
│   └── 04_integration/
│       └── prepare_magma_gene_sets.R
└── results/
    └── README.md
```

## Software
Gene and gene-set analyses were performed using MAGMA v1.10.
S-LDSC analyses were performed using the CBIIT implementation of LDSC (commit `1f09cf0`).
Detailed software dependencies are provided in `environment/environment.yml`.

## References

### Methods
- MAGMA: [de Leeuw et al. (2015), *PLoS Computational Biology*](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004219)
- LD Score regression: [Bulik-Sullivan et al. (2015), *Nature Genetics*](https://www.nature.com/articles/ng.3211)
- [Stratified LD Score regression: Finucane et al. (2015), *Nature Genetics*](https://www.nature.com/articles/ng.3404)

### GWAS
- ADHD - [Demontis et al. (2023), *Nature Genetics*](https://www.nature.com/articles/s41588-022-01285-8)
- Anxiety - [Storm et al. (2026), *Nature Genetics*](https://www.nature.com/articles/s41588-025-02485-8)
- ASD - [Grove et al. (2019), *Nature Genetics*](https://www.nature.com/articles/s41588-019-0344-8)
- BD - [O'Connell et al. (2025), *Nature*](https://www.nature.com/articles/s41586-024-08468-9)
- MDD - [Adams et al. (2025), *Cell*](https://www.cell.com/cell/fulltext/S0092-8674(24)01415-6)
- PTSD - [Nievergelt et al. (2024), *Nature Genetics*](https://www.nature.com/articles/s41588-024-01707-9)
- SCZ - [Trubetskoy et al. (2022), *Nature*](https://www.nature.com/articles/s41586-022-04434-5)

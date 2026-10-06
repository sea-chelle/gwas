# GWAS integration with transcriptomic and epigenomic data

This repository contains analysis workflows for integrating psychiatric genome-wide association study (GWAS) summary statistics with transcriptomic and epigenomic datasets using MAGMA and stratified LD score regression (S-LDSC). The pipeline supports gene-set and genomic annotation-based analyses for investigating relationships between GWAS genetic associations and features derived from RNA-seq and ATAC-seq data.

## Analyses

### MAGMA
MAGMA is used to perform gene-level association analysis and test enrichment of GWAS association within user-defined gene sets. Gene sets can be derived from transcriptomic data or regulatory regions linked to putative target genes.

### S-LDSC
S-LDSC is used to test whether genomic annotations derived from transcriptomic or epigenomic analyses explain disproportionate SNP heritability. Custom annotations are analysed conditional on the baselineLD v2.3 model. The workflow also supports joint S-LDSC models for directly comparing related annotations using regression coefficient covariance. Examples include sex (male vs female comparisons).

### Locus-level integration
Published fine-mapping results were integrated with transcriptomic and regulatory annotations from RNA-seq and ATAC-seq data to identify candidate psychiatric risk variants overlapping associated genomic features. Fine-mapped variants were defined according to the credible-set or candidate variant criteria reported by each original GWAS study.

Fine-mapped variants were intersected with:
- ±100 kb windows surrounding stress-responsive RNA-seq genes
- lifted human (GRCh37/hg19) coordinates of stress-responsive ATAC-seq regions

Intersections were performed separately for each annotation using BEDTools. 

| Trait | Study | Fine-mapping method / inclusion criterion | Variants |
|---|---|---|---:|
| ASD | Grove et al. 2019 | Variants in published CAVIAR credible sets | 380 |
| ADHD | Demontis et al. 2023 | Variants included in credible sets identified by all three methods: PAINTOR, CAVIARBF and FINEMAP | 1,139 |
| ANX | Strom et al. 2026 | Variants from six reported FINEMAP credible sets with configuration posterior probability >0.95 | 30 |
| MDD | Adams et al. 2025 | Variants in published PolyFun/SuSiE 95% credible causal sets | 14,652 |
| PTSD | Nievergelt et al. 2024 | Variants in published fine-mapping credible sets | 2,069 |
| SCZ | Trubetskoy et al. 2022 | Variants in FINEMAP 95% credible sets for 249 regions predicted to contain ≤3 causal variants (`k < 3.5`) | 20,591 |
| BD | O'Connell et al. 2025 | Published fine-mapped variants with PIP >0.5; full 95% credible-set membership was not publicly available | 80 |

Fine-mapping approaches and reporting criteria differed between studies. Variants were therefore selected according to the published credible-set or fine-mapped candidate definitions for each study rather than by imposing a uniform variant-level posterior probability threshold. Variant counts should consequently not be compared directly between disorders.

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
- S-LDSC: [Stratified LD Score regression: Finucane et al. (2015), *Nature Genetics*](https://www.nature.com/articles/ng.3404)

### GWAS
- ADHD: [Demontis et al. (2023), *Nature Genetics*](https://www.nature.com/articles/s41588-022-01285-8)
- ANX: [Storm et al. (2026), *Nature Genetics*](https://www.nature.com/articles/s41588-025-02485-8)
- ASD: [Grove et al. (2019), *Nature Genetics*](https://www.nature.com/articles/s41588-019-0344-8)
- BD: [O'Connell et al. (2025), *Nature*](https://www.nature.com/articles/s41586-024-08468-9)
- MDD: [Adams et al. (2025), *Cell*](https://www.cell.com/cell/fulltext/S0092-8674(24)01415-6)
- PTSD: [Nievergelt et al. (2024), *Nature Genetics*](https://www.nature.com/articles/s41588-024-01707-9)
- SCZ: [Trubetskoy et al. (2022), *Nature*](https://www.nature.com/articles/s41586-022-04434-5)

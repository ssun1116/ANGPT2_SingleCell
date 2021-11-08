# Transcriptomic feature of Ang2 (ANGPT2) contribution to T-cell exclusion in melanoma

Suggestion from Dr. Ben Izar
1. Determine the major cell type expressing angpt2
2. Develop a single cell signature associated with angpt2 high vs low 
3. Correlate the entire signature with bulk TCGA data to infer the associated T cell infiltration

## Project 1. Analysis of human scRNA-seq data (published)

1) Develop a single cell signature associated with angpt2 high vs. low from melanoma human scRNA-seq data 
 (n=31 including untreated, post-immunotherapy resistant, and post-immunotherapy responder)
* Code and expected results/output can be found here: https://github.com/livnatje/ImmuneResistance 

2) Correlate the entire signature with bulk TCGA data (TCGA-SKCM, metastasis vs primary in Sample Type) to
infer the association with T-cell tumor infiltration

## Project 2. Analysis of TCGA-SKCM data set

1) Stratify Ang2 high and low (e.g. >90th (high) vs. <50th (low) expression)2 using TCGA-SKCM data 
(Metastasis samples and Primary tumor samples, respectively)

2) Compare T-cell exclusion signatures1 in Ang2 high and Ang2 low groups from the TCGA-SKCM data. 

GSE115978 파일 구성
- GSE115978_cell.annotations.csv
- GSE115978_counts.csv
- GSE115978_tpm.csv : Transcripts per million.
 (step 1 : normalize for gene length, step 2 : normlize for sequencing depth)

**Paper link : https://doi.org/10.1016/j.cell.2018.09.006**

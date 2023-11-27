# Test of Polygenic Selection using Outlier Estimation

The given script uses a set of disease associated GWAS SNPs to estimate if the disease as a whole has undergone background slection (B-statistics) or recent positive selection when compared to control sets.

## Input Data
1. GWAS input file (chr  pos  allele  Beta  P-value)
2. Matched control SNPS (Input_SNP  control_1  control_2  ...  control_n)
3. Pre computed B - statistics (McVicker et al., 2009)
4. Pre computed integrated Haplotype Score (iHS) (Johnson et al., 2018)

## Ouput file
1. bstat_result.txt - Reports GWAS AUC percentile when comapred to controls
2. iHS_Pop_result.txt - Reports GWAS AUC percentile when comapred to controls for the specific 1000 genomes population

#!/bin/bash

# curated clinical data
wget https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp.gz
gzip -d Survival_SupplementalTable_S1_20171025_xena_sp.gz

# immune subtype
wget https://pancanatlas.xenahubs.net/download/Subtype_Immune_Model_Based.txt.gz
gzip -d Subtype_Immune_Model_Based.txt.gz

# molecular subtype
wget https://pancanatlas.xenahubs.net/download/TCGASubtype.20170308.tsv.gz
gzip -d TCGASubtype.20170308.tsv.gz

# sample type and primary disease
wget https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
gzip -d TCGA_phenotype_denseDataOnlyDownload.tsv.gz


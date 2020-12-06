#!/bin/bash

# RSEM expected count, log2(x+1) transformed
wget https://toil.xenahubs.net/download/tcga_gene_expected_count.gz
gzip -d tcga_gene_expected_count.gz

# TPM, log2(x+0.001) transformed
wget https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz
gzip -d tcga_RSEM_gene_tpm.gz

# gene info
wget https://toil.xenahubs.net/download/probeMap/gencode.v23.annotation.gene.probemap


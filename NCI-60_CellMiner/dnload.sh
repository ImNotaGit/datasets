#!/bin/bash

# cell line metadata
wget https://discover.nci.nih.gov/cellminer/samples/NCI60_CELL_LINE_METADATA.txt -P metadata

# exome mutations
wget https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_DNA__Exome_Seq_none.zip
unzip nci60_DNA__Exome_Seq_none.zip
mv output exome

# mRNA (5 platforms merged)
wget https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__5_Platform_Gene_Transcript_Average_z_scores.zip
unzip nci60_RNA__5_Platform_Gene_Transcript_Average_z_scores.zip
mv output rna_ave_zscore

# mRNA (RNA-seq gene level)
wget https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__RNA_seq_composite_expression.zip
unzip nci60_RNA__RNA_seq_composite_expression.zip
mv output rnaseq

# mRNA (RNA-seq isoform level)
wget https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__RNA_seq_isoforms.zip
unzip nci60_RNA__RNA_seq_isoforms.zip
mv output rnaseq_isoform

# miRNA
wget https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__microRNA_OSU_V3_chip_log2.zip
unzip nci60_RNA__microRNA_OSU_V3_chip_log2.zip
mv output mirna

# protein array
wget https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_Protein__Lysate_Array_log2.zip
unzip nci60_Protein__Lysate_Array_log2.zip
mv output protein_array

# compound activity
# raw
wget https://discover.nci.nih.gov/cellminerdata/rawdata/DTP_NCI60_RAW.zip
unzip DTP_NCI60_RAW.zip
mv output compound_raw
# processed
wget https://discover.nci.nih.gov/cellminerdata/normalizedArchives/DTP_NCI60_ZSCORE.zip
unzip DTP_NCI60_ZSCORE.zip
mv output compound_zscore

# metabolomics
wget https://wiki.nci.nih.gov/download/attachments/155845004/WEB_DATA_METABOLON.ZIP
unzip WEB_DATA_METABOLON.ZIP -d metabolome


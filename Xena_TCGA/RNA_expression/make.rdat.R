library(my.utils)

##### 1. based on RSEM expected count

dat <- fread("tcga_gene_expected_count")
gid <- dat$sample
mat <- data.matrix(dat[,-1])
colnames(mat) <- tolower(colnames(mat))
rownames(mat) <- gid
mat <- 2^mat-1 # downloaded values are log2(x+1) transformed
mat[mat<0] <- 0 # make sure no negative values due to numerical error

gene.info <- fread("gencode.v23.annotation.gene.probemap")
setnames(gene.info, c("id","symbol","chr","start","end","strand"))
nrow(mat)==nrow(gene.info) # TRUE
gene.info <- gene.info[match(rownames(mat), id)]
all(rownames(mat)==gene.info$id) # TRUE

phe <- readRDS("../phenotype/pheno.RDS")
all(colnames(mat) %in% phe$sample.id) # FALSE
sum(!colnames(mat) %in% phe$sample.id) # 1
phe <- phe[match(colnames(mat), sample.id)]
idx <- phe[, which(is.na(sample.id))]
phe[idx, sample.id:=colnames(mat)[idx]]
phe[idx, patient.id:=str_sub(sample.id, 1, 12)]
phe <- phe[order(cancer.type)]
mat <- mat[, match(phe$sample.id, colnames(mat))]
all(phe$sample.id==colnames(mat)) # TRUE

# RSEM raw count (NOT on the log scale), all samples (including all sample types, e.g. tumor-adjacent normal tissues)
tcga.xena.rsem.count <- list(expr=mat, pheno=phe, geneid=gene.info)
saveRDS(tcga.xena.rsem.count, file="tcga.xena.all.rsem.count.RDS")

# TMM normalize across all samples (i.e. all cancer types), including all sample types e.g. tumor-adjacent normal tissues
res <- rm.low.genes(tcga.xena.rsem.count, rm.low.frac.gt=0.995)
res$expr <- get.tmm.log.cpm(res$expr)
saveRDS(res, file="tcga.xena.all.tmm.log.cpm.RDS")

# TMM normalize across all samples (i.e. all cancer types) with sample type 1 or 3 (primary solid tumor or blood-derived cancer samples)
res <- filter.eset(tcga.xena.rsem.count, sample.type.id %in% c(1,3))
res <- rm.low.genes(res, rm.low.frac.gt=0.995)
res$expr <- get.tmm.log.cpm(res$expr)
saveRDS(res, file="tcga.xena.primary.tumor.tmm.log.cpm.RDS")

# TMM normalize by cancer type, each cancer type includes all sample types e.g. tumor-adjacent normal tissues
cts <- unique(res$pheno$cancer.type)
names(cts) <- cts
cts <- cts[!is.na(cts)]
res <- lapply(cts, function(x) {
  tmp <- filter.eset(tcga.xena.rsem.count, cancer.type==x)
  tmp <- rm.low.genes(tmp, rm.low.frac.gt=0.9)
  tmp$expr <- get.tmm.log.cpm(tmp$expr)
  tmp
})
saveRDS(res, file="tcga.xena.all.by.cancer.type.tmm.log.cpm.RDS")

# TMM normalize by cancer type, each cancer type includes only sample type 1 or 3 (primary solid tumor or blood-derived cancer samples)
cts <- unique(res$pheno$cancer.type)
names(cts) <- cts
cts <- cts[!is.na(cts)]
res <- lapply(cts, function(x) {
  tmp <- filter.eset(tcga.xena.rsem.count, cancer.type==x & sample.type.id %in% c(1,3))
  tmp <- rm.low.genes(tmp, rm.low.frac.gt=0.9)
  tmp$expr <- get.tmm.log.cpm(tmp$expr)
  tmp
})
saveRDS(res, file="tcga.xena.primary.tumor.by.cancer.type.tmm.log.cpm.RDS")


##### 2. based on TPM

dat <- fread("tcga_RSEM_gene_tpm")
gid <- dat$sample
mat <- data.matrix(dat[,-1])
colnames(mat) <- tolower(colnames(mat))
rownames(mat) <- gid
mat <- 2^mat-0.001 # downloaded values are log2(x+0.001) transformed
mat[mat<0] <- 0 # make sure no negative values due to numerical error

gene.info <- fread("gencode.v23.annotation.gene.probemap")
setnames(gene.info, c("id","symbol","chr","start","end","strand"))
nrow(mat)==nrow(gene.info) # TRUE
gene.info <- gene.info[match(rownames(mat), id)]
all(rownames(mat)==gene.info$id) # TRUE

phe <- readRDS("../phenotype/pheno.RDS")
all(colnames(mat) %in% phe$sample.id) # FALSE
sum(!colnames(mat) %in% phe$sample.id) # 1
phe <- phe[match(colnames(mat), sample.id)]
idx <- phe[, which(is.na(sample.id))]
phe[idx, sample.id:=colnames(mat)[idx]]
phe[idx, patient.id:=str_sub(sample.id, 1, 12)]
phe <- phe[order(cancer.type)]
mat <- mat[, match(phe$sample.id, colnames(mat))]
all(phe$sample.id==colnames(mat)) # TRUE

# TPM values (NOT on the log scale), all samples
tcga.xena.tpm <- list(expr=mat, pheno=phe, geneid=gene.info)
saveRDS(tcga.xena.tpm, file="tcga.xena.all.tpm.RDS")



library(my.utils)


##### cell line meta-data

phe <- fread("metadata/data.tsv")
phe <- phe[,1:12]
setnames(phe, c("name","tissue","age","gender","prior.trt","is.epithelial","histology","source","ploidy","p53","mdr","doubling.time"))


##### mRNA expression, 5 platform transcript average z-score

dat <- fread("rna_ave_zscore/data.tsv")

genes <- dat[,1:9]
setnames(genes, c("probe.id","symbol","entrez.id","chr","start","end","cytoband","mrna.refseq","prot.refseq"))

dat <- dat[,-1:-9]
names(dat) <- phe$name # I've checked that the orders are the same
dat <- data.matrix(dat)
rownames(dat) <- genes$symbol

nci60.exprs.z <- list(expr=dat, pheno=phe, geneid=genes)
saveRDS(nci60.exprs.z, file="rna.ave.zscore.RDS")


##### mRNA expression: RNA-seq, log2(FPKM+1)

dat <- fread("rnaseq/data.tsv")

genes <- dat[,1:6]
setnames(genes, c("symbol","entrez.id","chr","start","end","cytoband"))

dat <- dat[,-1:-6]
names(dat) <- phe$name # I've checked that the orders are the same - although the spelling/format is slightly different
dat <- data.matrix(dat)
rownames(dat) <- genes$symbol

# values are log2(FPKM+1), transform into log2(TPM+1)
dat <- 2^(dat)-1
dat <- apply(dat, 2, function(x) x/sum(x)*1e6)
dat <- log2(dat+1)

nci60.log.tpm <- list(expr=dat, pheno=phe, geneid=genes)
saveRDS(nci60.log.tpm, file="rna.log.tpm.RDS")


##### metabolome, averaged value of 3 replicates

tmp <- fread("metabolome/WEB_DATA_METABOLON.TXT")
mets <- unique(tmp[,1:2])
setnames(mets, c("id","name"))
mat <- dcast(tmp, MOLTID ~ cellname, value.var="VALUE")
tmp <- mat$MOLTID
mat <- data.matrix(mat[,-1])
rownames(mat) <- tmp
all(rownames(mat)==mets$id) # TRUE

tmp1 <- str_replace_all(colnames(mat), "-", "_")
tmp1 <- str_replace_all(tmp1, " ", "")
tmp2 <- str_split(phe$name, ":", simplify=TRUE)[,2]
# checking:
setdiff(tmp1, tmp2)
setdiff(tmp2, tmp1)
# fix cell line names:
tmp1[tmp1=="A549/ATCC"] <- "A549"
tmp1[tmp1=="HL_60(TB)"] <- "HL_60"
tmp1[tmp1=="MDA_MB_231/ATCC"] <- "MDA_MB_231"
tmp1[tmp1=="NCI/ADR_RES"] <- "NCI_ADR_RES"
tmp1[tmp1=="T_47D"] <- "T47D"
# now:
all(tmp1 %in% tmp2) # TRUE
colnames(mat) <- phe$name[match(tmp1, tmp2)]
phe <- phe[name %in% colnames(mat)]
mat <- mat[, phe$name]
nci60.mets <- list(expr=mat, pheno=phe, geneid=mets)
saveRDS(nci60.mets, file="metabolome.RDS")



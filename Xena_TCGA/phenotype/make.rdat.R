library(data.table)
library(stringr)

d1 <- fread("Survival_SupplementalTable_S1_20171025_xena_sp") # for details see https://doi.org/10.1016/j.cell.2018.02.052
d1 <- d1[, lapply(.SD, function(x) if (is.character(x)) tolower(x) else x)]
d1 <- d1[, c(1:6,26:33,7:10,15,18:21,23,24,25,12,34)]
setnames(d1, c("sample.id","patient.id","cancer.type","age","gender","race","os","os.days","dss","dss.days","dfi","dfi.days","pfi","pfi.days","patho.stage","clin.stage","histo.type","histo.grade","tumor.stat","cod","new.tumor.type","new.tumor.site","new.tumor.site1","mor.first.course","margin.stat","resid.tumor","menopause","redaction"))
d1 <- d1[, lapply(.SD, function(x) ifelse(x %in% c("","[discrepancy]","[unknown]","[not evaluated]"), NA, x))]
d1[startsWith(patho.stage, "stage"), patho.stage:=str_sub(patho.stage, 7)]
d1[startsWith(clin.stage, "stage"), clin.stage:=str_sub(clin.stage, 7)]
d1[histo.grade=="high grade", histo.grade:="high"]
d1[histo.grade=="low grade", histo.grade:="low"]
d1[mor.first.course=="complete remission/response", mor.first.course:="cr"]
d1[mor.first.course=="partial remission/response", mor.first.course:="pr"]
d1[mor.first.course=="progressive disease", mor.first.course:="pd"]
d1[mor.first.course=="stable disease", mor.first.course:="sd"]

d2 <- fread("TCGA_phenotype_denseDataOnlyDownload.tsv")
d2 <- d2[, lapply(.SD, function(x) if (is.character(x)) tolower(x) else x)]
setnames(d2, c("sample.id","sample.type.id","sample.type","primary.disease"))
d2[sample.type=="", sample.type:=NA]
phe <- merge(d2, d1, by="sample.id", all=TRUE)

d3 <- fread("TCGASubtype.20170308.tsv")
d3 <- d3[, lapply(.SD, function(x) if (is.character(x)) tolower(x) else x)]
setnames(d3, str_replace_all(tolower(names(d3)), "_", "."))
names(d3)[c(1,3)] <- c("sample.id","subtype.dna.meth")
d3 <- d3[, lapply(.SD, function(x) ifelse(x %in% c("","-","#n/a"), NA, x))]
phe <- merge(phe, d3, by="sample.id", all=TRUE)

d4 <- fread("Subtype_Immune_Model_Based.txt")
d4 <- d4[, lapply(.SD, function(x) if (is.character(x)) tolower(x) else x)]
setnames(d4, c("sample.id","subtype.immune"))
phe <- merge(phe, d4, by="sample.id", all=TRUE)

phe[is.na(patient.id), patient.id:=str_sub(sample.id, 1, 12)]
mapp <- c("uterine corpus endometrioid carcinoma"="ucec", "kidney papillary cell carcinoma"="kirp", "thymoma"="thym", "stomach adenocarcinoma"="stad", "kidney clear cell carcinoma"="kirc", "adrenocortical cancer"="acc", "glioblastoma multiforme"="gbm", "ovarian serous cystadenocarcinoma"="ov", "lung squamous cell carcinoma"="lusc", "bladder urothelial carcinoma"="blca", "esophageal carcinoma"="esca", "liver hepatocellular carcinoma"="lihc", "sarcoma"="sarc", "mesothelioma"="meso", "skin cutaneous melanoma"="skcm", "thyroid carcinoma"="thca", "acute myeloid leukemia"="laml", "diffuse large b-cell lymphoma"="dlbc", "kidney chromophobe"="kich", "pheochromocytoma & paraganglioma"="pcpg", "uveal melanoma"="uvm", "lung adenocarcinoma"="luad", "prostate adenocarcinoma"="prad", "testicular germ cell tumor"="tgct", "pancreatic adenocarcinoma"="paad", "cervical & endocervical cancer"="cesc", "breast invasive carcinoma"="brca", "colon adenocarcinoma"="coad", "cholangiocarcinoma"="chol", "head & neck squamous cell carcinoma"="hnsc", "rectum adenocarcinoma"="read", "brain lower grade glioma"="lgg", "uterine carcinosarcoma"="ucs")
phe[is.na(cancer.type), cancer.type:=mapp[primary.disease]]
phe <- phe[order(cancer.type)]

saveRDS(phe, file="pheno.RDS")


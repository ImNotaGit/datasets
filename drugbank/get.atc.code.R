library(xml2)
library(parallel)
library(data.table)
library(stringr)
nc <- 32L

# get the mapping between drugs and their ATC codes, create drug sets based on the different levels of ATC code

# read in xml
dbxml <- read_xml("drugbank.xml")
# the number of children nodes, i.e. the number of drugs
xml_length(dbxml) # 13475, as of drugbank release version 5.1.5

# as the xml file has a namespace, each tag name need to be prefixed by "d1:", representing the default name space.
# if the root node is simply "<drugbank>" w/o namespace definition, then drop the "d1:"
# drugbank primary id
id <- xml_text(xml_find_all(dbxml, "./d1:drug/d1:drugbank-id[@primary='true']"))
# drug name
name <- xml_text(xml_find_all(dbxml, "./d1:drug/d1:name"))
# ATC codes
ii <- 1:xml_length(dbxml)
xx <- xml_children(dbxml)
atc <- mclapply(ii, function(i) {
  x <- xx[[i]]
  message(i)
  did <- xml_attr(xml_find_all(x, "./d1:atc-codes/d1:atc-code"), "code") # drug-level ATC codes
  ids <- xml_attr(xml_find_all(x, "./d1:atc-codes/d1:atc-code/d1:level"), "code") # higher level ATC codes
  tms <- xml_text(xml_find_all(x, "./d1:atc-codes/d1:atc-code/d1:level")) # higher level terms
  list(did, data.table(atc.code=ids, atc.term=tms))
}, mc.cores=nc)

# create drug-level ATC code mapping (each drug can have multiple drug-level ATC codes)
drug.level.atc <- data.table(drug.id=id, drug.name=name,
                             atc.code=lapply(atc, function(x) {
                                        if (length(x[[1]])==0) res <- NA else res <- x[[1]]
                                        res
                                      })
)

names(atc) <- id
atc1 <- rbindlist(lapply(atc, function(x) x[[2]]), idcol="drug.id")
atc1[, drug.name:=name[match(drug.id, id)]]

# create a mapping between ATC code and term for higher-level ATC codes
atc.anno <- unique(atc1[, .(atc.code, atc.term)])[order(atc.code)]

# create drug set based on ATC
x <- atc.anno$atc.code # use ATC code because the terms are not unique (multiple codes can have the same term)
names(x) <- x
# using drug name
atc.drug.set <- lapply(x, function(i) atc1[atc.code==i, unique(drug.name)])
# using drug id
atc.drug.id.set <- lapply(x, function(i) atc1[atc.code==i, unique(drug.id)])


save(atc.anno, drug.level.atc, atc.drug.set, atc.drug.id.set, file="drug.atc.RData")




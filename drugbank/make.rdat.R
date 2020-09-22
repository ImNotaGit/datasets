library(xml2)
library(parallel)
library(data.table)
library(stringr)
nc <- 32L

# I decided to save the drug info and drug-protein mapping data separately, parsed from the full database (.xml file). I will clean the protein info from the downloaded protein info files (.csv files).

# read in xml
dbxml <- read_xml("drugbank.xml")
# the number of children nodes, i.e. the number of drugs
xml_length(dbxml) # 13580, as of drugbank release version 5.1.7


## ----drugs info----

# as the xml file has a namespace, each tag name need to be prefixed by "d1:", representing the default name space.
# if the root node is simply "<drugbank>" w/o namespace definition, then drop the "d1:"
# drugbank primary id
id <- xml_text(xml_find_all(dbxml, "./d1:drug/d1:drugbank-id[@primary='true']"))
# drug type: "biotech" or "small molecule"
type <- xml_attr(xml_children(dbxml), "type")
# groups: e.g. "approved", "experimental", etc. If belongs to multiple groups, the groups are joint by ; into a single string
group <- sapply(xml_find_all(dbxml, "./d1:drug/d1:groups"), function(x) paste(xml_text(xml_children(x)), collapse=" ; "))
# drug name
name <- xml_text(xml_find_all(dbxml, "./d1:drug/d1:name"))
# CAS number
cas <- xml_text(xml_find_all(dbxml, "./d1:drug/d1:cas-number"))
# drug synonyms
synonym <- sapply(xml_find_all(dbxml, "./d1:drug/d1:synonyms"), function(x) paste(unique(xml_text(xml_children(x))), collapse=" ; "))
# drug product names, including international brand names. This kind of operation (xml_find_all within apply) is very slow so I do it in parallel:
product <- unlist(mclapply(xml_children(dbxml), function(x) paste(unique(c(xml_text(xml_find_all(x, "./d1:products/d1:product/d1:name")), xml_text(xml_find_all(x, "./d1:international-brands/d1:international-brand/d1:name")))), collapse=" ; "), mc.cores=nc))
# drug product labeller (drug company)
labeller <- unlist(mclapply(xml_children(dbxml), function(x) paste(unique(xml_text(xml_find_all(x, "./d1:products/d1:product/d1:labeller"))), collapse=" ; "), mc.cores=nc))
# indication
indication <- xml_text(xml_find_all(dbxml, "./d1:drug/d1:indication"))
# mechanism of action
moa <- xml_text(xml_find_all(dbxml, "./d1:drug/d1:mechanism-of-action"))

# save drugs info
drugs <- data.table(id=id, type=type, group=group, cas=cas, name=name, synonym=synonym, product=product, labeller=labeller, indication=indication, moa=moa)
drugs[drugs==""] <- NA
saveRDS(drugs, file="drugs.RDS")


## ----drug-protein mapping----

# the proteins include targets, enzymes, carriers and transporters. For now, get the targets only
mapdrug <- function(drug) {
  tryCatch(
    {
      # drug name
      drugname <- tolower(xml_text(xml_find_all(drug, "./d1:name")))
      # CAS number
      cas <- xml_text(xml_find_all(drug, "./d1:cas-number"))
      # drugbank primary id
      drugid <- xml_text(xml_find_all(drug, "./d1:drugbank-id[@primary='true']"))
      # targets
      target <- lapply(xml_find_all(drug, "./d1:targets/d1:target"), function(x) {
        id <- xml_text(xml_find_all(x, "./d1:id"))
        if (length(id)==0) id <- NA else if (id=="") id <- NA
        uniprot <- xml_attr(xml_find_all(x, "./d1:polypeptide"), "id")
        if (length(uniprot)==0) uniprot <- NA else if (uniprot=="") uniprot <- NA
        symbol <- xml_text(xml_find_all(x, "./d1:polypeptide/d1:gene-name"))
        if (length(symbol)==0) symbol <- NA else if (symbol=="") symbol <- NA
        action <- paste(unique(tolower(xml_text(xml_find_all(x, "./d1:actions/d1:action")))), collapse=" ; ")
        if (length(action)==0) action <- NA else if (action=="") action <- NA
        pharmaco.active <- xml_text(xml_find_all(x, "./d1:known-action"))
        if (length(pharmaco.active)==0) pharmaco.active <- NA else if (pharmaco.active=="") pharmaco.active <- NA
        data.table(protein.id=id, protein.uniprot.id=uniprot, protein.gene.symbol=symbol, action=action, pharmaco.active=pharmaco.active)
      })
      if (length(target)==0) {
        target <- data.table(protein.id=NA, protein.uniprot.id=NA, protein.gene.symbol=NA, action=NA, pharmaco.active=NA)
      } else target <- do.call(rbind, target)
      cbind(data.table(drug.id=drugid, drug.name=drugname, cas=cas, protein.type="target"), target)
    },
    error=function(e) return(e)
  )
}
tmp <- mclapply(xml_children(dbxml), mapdrug, mc.cores=nc)
save(tmp, file="tmp.RData")
# rbind and save tmp if needed
tryCatch(
  drugs2targets <- do.call(rbind, tmp),
  error=function(e) {
    print(e)
    save(tmp, file="tmp.RData")
 },
  warning=function(w) {
    print(w)
    save(tmp, file="tmp.RData")
  }
)
# a summary of actions
action.simp <- drugs2targets[, .N, by=action][order(-N)]
write.table(action.simp, file="interpret_actions.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# **manually edit** interpret_actions.tsv to add simplified action (activation/inhibition) for each action term, then re-readin
action.simp <- fread("interpret_actions.tsv")
# make two vectors as dictionaries
act.strict <- action.simp[, action.simp.strict]
names(act.strict) <- action.simp[, action]
act.loose <- action.simp[, action.simp.loose]
names(act.loose) <- action.simp[, action]
save(action.simp, act.strict, act.loose, file="interpret_actions.RData")
# add simplified action to drugs2targets
drugs2targets[, c("action.simp.strict", "action.simp.loose"):=list(act.strict[action], act.loose[action])]
setcolorder(drugs2targets, c(names(drugs2targets)[1:8], "action.simp.strict", "action.simp.loose", "pharmaco.active"))
saveRDS(drugs2targets, file="drugs2targets.RDS")


## ----protein info----

# targets
targets <- fread("targets/all.csv")
setnames(targets, str_replace_all(tolower(names(targets)), " ", "\\."))
setnames(targets, "gene.name", "gene.symbol")
targets <- as.data.table(lapply(targets, function(x) str_replace_all(x, "; ", " ; ")))
saveRDS(targets, file="targets.RDS")



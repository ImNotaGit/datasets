library(data.table)
library(stringr)


## collect gene sets

make.gsets <- function(fn, name=1) {
  gsets <- readLines(fn)
  gsets <- str_split(gsets, "\t")
  names(gsets) <- sapply(gsets, function(x) x[name]) # gene set name in the 1st position; 2nd position is ID
  lapply(gsets, function(x) x[-1:-2]) # gene list starts from the 3rd position; change as appropriate
}
reactome <- make.gsets("ReactomePathways.gmt")
saveRDS(reactome, file="reactome.gsets.RDS")

# also name by ID (used for subsetting later, as ID may be more reliable than name; r1 is in the same order as reactome)
r1 <- make.gsets("ReactomePathways.gmt", 2)


## collect subsets of immune-related and metabolism-related gene sets

# pathway list
paths <- fread("ReactomePathways.txt", header=FALSE)
setnames(paths, c("id","name","species"))
paths <- paths[species=="Homo sapiens"]
# pathway hierarchical structure
hier <- fread("ReactomePathwaysRelation.txt", header=FALSE)
setnames(hier, c("parent","child"))
# function to get all descendants of given node(s), including the given node(s)
get.desc <- function(top.terms) {
  top.ids <- paths[name %in% top.terms, id]
  res <- top.ids
  cur.ids <- hier[parent %in% top.ids, child]
  while (length(cur.ids)!=0) {
    res <- c(res, cur.ids)
    cur.ids <- hier[parent %in% cur.ids, child]
  }
  res
}

# immune-related
immune.ids <- get.desc("Immune System")
immune.gsets <- reactome[names(r1) %in% immune.ids]
saveRDS(immune.gsets, file="reactome.immune.gsets.RDS")
# metabolism-related
# 1.
metab.ids <- get.desc(c("Metabolism","Transport of small molecules"))
metab.gsets <- reactome[names(r1) %in% metab.ids]
saveRDS(metab.gsets, file="reactome.metabolism.gsets.RDS")
# 2. with "Metabolism of proteins","Metabolism of RNA"
metab.ids <- get.desc(c("Metabolism","Metabolism of proteins","Metabolism of RNA","Transport of small molecules"))
metab.gsets <- reactome[names(r1) %in% metab.ids]
saveRDS(metab.gsets, file="reactome.metabolism.gsets.extra.RDS")


library(stringr)

make.gsets <- function(fn, name=1) {
  gsets <- readLines(fn)
  gsets <- str_split(gsets, "\t")
  names(gsets) <- sapply(gsets, function(x) x[name]) # gene set name in the 1st position
  lapply(gsets, function(x) x[-1:-2]) # gene list starts from the 3rd position; change as appropriate
}

fns <- dir(pattern=".*\\.gmt")
for (fn in fns) {
  res <- make.gsets(fn)
  fo <- paste0(str_sub(fn, 1, -5), ".RDS")
  saveRDS(res, file=fo)
}

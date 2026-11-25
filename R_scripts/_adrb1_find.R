suppressMessages(library(TwoSampleMR))
ao <- try(available_outcomes(), silent = TRUE)
if (inherits(ao, "try-error")) { cat("ERR:", conditionMessage(attr(ao,"condition")), "\n"); quit(status=1) }
saveRDS(ao, "ao_all.rds")
cat("total:", nrow(ao), "\n")
g <- ao[grepl("ENSG00000043591", ao$id) | grepl("\\bADRB1\\b", ao$trait, ignore.case=TRUE), ]
cat("hits:", nrow(g), "\n")
if (nrow(g)) print(g[order(-g$sample_size), c("id","trait","sample_size","nsnp")][1:min(10,nrow(g)),], row.names=FALSE)

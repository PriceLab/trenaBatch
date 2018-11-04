# creates aggregate of all models in new data.frametbl.models
# adding hugo geneSymbols for each tf, ordering by abs(pearsonCoeff)
#-------------------------------------------------------------------------------------------------------
target.dir <- "~/github/trenaBatch/placenta/footprints/out"
all.results <- list.files(target.dir)
all.results <- grep(".out.RData", all.results, value=TRUE)
printf("gene result count: %d", length(all.results)) # 14099

tbls.all <- list()
#max <- 3 # length(all.gene.directories)
max <- length(all.results)

for(filename in all.results[1:max]){
   targetGene <- sub(".out.RData", "", filename, fixed=TRUE)
   full.path <- file.path(target.dir, filename)
   stopifnot(file.exists(full.path))
   load(full.path)
   tbl <- x$model
   if(nrow(tbl) == 0) next
   tbl <- tbl[order(abs(tbl$pearsonCoeff), decreasing=TRUE),]       
   if(nrow(tbl) > 20)
      tbl <- tbl[1:20,]
   tbl$targetGene <- targetGene
   tbl$rank <- seq_len(nrow(tbl))
   printf("model for %s: %d rows.orig, %d rows.trimmed, %d cols", targetGene, nrow(x$model), nrow(tbl), ncol(tbl))
   tbls.all[[targetGene]] <- tbl
   } # for filename

tbl.models <- do.call(rbind, tbls.all)
rownames(tbl.models) <- NULL

save(tbl.models, file="tbl.models.all.placenta.footprints.RData")
    

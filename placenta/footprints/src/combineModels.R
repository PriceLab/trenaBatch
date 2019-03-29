# creates aggregate of all models in new data.frametbl.models
# adding hugo geneSymbols for each tf, ordering by abs(pearsonCoeff)
#-------------------------------------------------------------------------------------------------------
filenames <- system("find 2019mar28/ -name models.RData", intern=TRUE)
printf("gene result count: %d", length(filenames)) # 14

tbls.all <- list()
max <- 3 # length(all.gene.directories)
#max <- length(all.results)

for(full.path in filenames){
   gene.gene <- strsplit(full.path, "/", fixed=TRUE)[[1]][2]
   targetGene <- strsplit(gene.gene, "-", fixed=TRUE)[[1]][1]
   stopifnot(file.exists(full.path))
   print(load(full.path))
   tbl <- tbls$model
   if(nrow(tbl) == 0) next
   tbl <- tbl[order(abs(tbl$pearsonCoeff), decreasing=TRUE),]       
   #if(nrow(tbl) > 20)
   #   tbl <- tbl[1:20,]
   tbl$targetGene <- targetGene
   tbl$rank <- seq_len(nrow(tbl))
   printf("model for %s: %d rows.orig, %d rows.trimmed, %d cols", targetGene, nrow(tbls$model), nrow(tbl), ncol(tbl))
   tbls.all[[targetGene]] <- tbl
   } # for full.path

tbl.models <- do.call(rbind, tbls.all)
rownames(tbl.models) <- NULL
colnames(tbl.models)[1] <- "tf"

save(tbl.models, file="tbl.models.demo.200.v2.placenta.old.matrix.footprints.RData")
    

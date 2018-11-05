# creates aggregate of all models in tbl.models, file is tbl.models.all.RData
# adding hugo geneSymbols for each tf, and for the target gene
#-------------------------------------------------------------------------------------------------------
target.dir <- "~/github/trenaSGM/inst/batch/cory/AD/out.2018-07-08.17:16:33"
all.gene.directories <- list.files(target.dir)
all.gene.directories <- grep("ENSG0", all.gene.directories, value=TRUE)
printf("gene directory count: %d", length(all.gene.directories))

tbls.all <- list()
max <- 3 # length(all.gene.directories)
max <- length(all.gene.directories)
for(gene.dir in all.gene.directories[1:max]){
    path.models <- file.path(target.dir, gene.dir, "models.RData")
    path.footprints <- file.path(target.dir, gene.dir, "tbl.fp.RData")
      # existence of footprints indicates usable model - i think
    if(file.exists(path.footprints)){
        load(path.models)
        tbl <- tbls$model
        if(nrow(tbl) == 0) next
        #print(gene.dir)
        if(gene.dir == "ENSG00000005001-PRSS22") browser()
        tbl <- tbl[order(tbl$pcaMax, decreasing=TRUE),]       
        tokens <- strsplit(gene.dir, "-")[[1]]
        ensembl.id <- tokens[1]
        if(length(tokens) == 2)
           genesym <- tokens[2]
         else
           genesym <- ensembl.id
        tbl$targetGene <- ensembl.id
        tbl$rank <- seq_len(nrow(tbl))
        printf("model for %s (%s): %d rows, %d cols", genesym, ensembl.id, nrow(tbls$model), ncol(tbl))
        tbls.all[[ensembl.id]] <- tbl
        } # if file.exists
     } # for gene.dir

  
tbl.models <- tbl <- do.call(rbind, tbls.all)
rownames(tbl.models) <- NULL

library(org.Hs.eg.db)
all.ensg <- sort(unique(c(tbl.models$gene, tbl.models$targetGene)))
tbl.map <- select(org.Hs.eg.db, all.ensg, "ENSEMBL", columns="SYMBOL")
dups <- which(duplicated(tbl.map$ENSEMBL))
if(length(dups) > 0)
    tbl.map <- tbl.map[-dups,]

NAs <- which(is.na(tbl.map$SYMBOL))
if(length(NAs) > 0)
    tbl.map$SYMBOL[NAs] <- tbl.map$ENSEMBL[NAs]

naStrings <- grep("NA", tbl.map$SYMBOL)
if(length(naStrings) > 0)
    tbl.map$SYMBOL[naStrings] <- tbl.map$ENSEMBL[naStrings]


map.list <- tbl.map$SYMBOL
names(map.list) <- tbl.map$ENSEMBL

tbl.models$tf.hgnc <- map.list[tbl.models$gene]
tbl.models$targetGene.hgnc <- map.list[tbl.models$targetGene]
colnames(tbl.models)[grep("^gene$", colnames(tbl.models))] <- "tf.ensembl"
colnames(tbl.models)[grep("^targetGene$", colnames(tbl.models))] <- "targetGene.ensembl"
coi <- c("targetGene.hgnc",
         "tf.hgnc",
         "targetGene.ensembl",
         "tf.ensembl",
         "betaLasso",
         "lassoPValue",
         "pearsonCoeff",
         "rfScore",
         "betaRidge",
         "spearmanCoeff",
         "betaSqrtLasso",
         "concordance",
         "pcaMax",
         "bindingSites",
         "rank")
tbl.models <- tbl.models[, coi]
save(tbl.models, file="tbl.models.all.RData")
    

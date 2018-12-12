# creates aggregate of all models in tbl.models, file is tbl.models.all.RData
# adding hugo geneSymbols for each tf, and for the target gene
#-------------------------------------------------------------------------------------------------------
target.dir <- "/ssd/cory/github/trenaBatch/AD/footprints/demo"
all.gene.directories <- list.files(target.dir)
all.gene.directories <- grep("ENSG0", all.gene.directories, value=TRUE)
printf("gene directory count: %d", length(all.gene.directories))

tbls.all <- list()
#max <- 10 # length(all.gene.directories)
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
        tbl <- tbl[order(tbl$rfScore, decreasing=TRUE),]       
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

library(TrenaProject)
load(system.file(package="TrenaProject", "extdata", "geneInfoTable.RData"))
dim(tbl.geneInfo)

colnames(tbl.models)[grep("^gene$", colnames(tbl.models))] <- "tf.ensg"
colnames(tbl.models)[grep("^targetGene$", colnames(tbl.models))] <- "target.ensg"

matches <- match(tbl.models$tf.ensg, tbl.geneInfo$ensg)
tbl.models$tf.symbol <- tbl.geneInfo$geneSymbol[matches]

matches <- match(tbl.models$target.ensg, tbl.geneInfo$ensg)
tbl.models$target.symbol <- tbl.geneInfo$geneSymbol[matches]
save(tbl.models, file="tbl.models.all.RData")
    

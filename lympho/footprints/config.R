library(TrenaProjectLymphocyte)

stopifnot(packageVersion("TrenaProject") >= "0.99.36")
stopifnot(packageVersion("TrenaProjectLymphocyte") >= "0.99.8")

OUTPUTDIR <- "2019feb27"
WORKERS <- 10

SOLVERS <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "sqrtLasso")

fp.logDir <- file.path(OUTPUTDIR, "logs.fp")
tfMapping.logDir <- file.path(OUTPUTDIR, "logs.tfMapping")
model.logDir <- file.path(OUTPUTDIR, "logs.model")

trenaProject <- TrenaProjectLymphocyte()
footprint.databases <- getFootprintDatabaseNames(trenaProject) # "lymphoblast_hint_20"
stopifnot(all(footprint.databases %in% getFootprintDatabaseNames(trenaProject)))

getExpressionMatrixNames(trenaProject)
matrix.name <- "GTEx.lymphocyte.ensg.matrix.asinh"
stopifnot(matrix.name %in% getExpressionMatrixNames(trenaProject))
mtx <- getExpressionMatrix(trenaProject, matrix.name)

print(load(system.file(package="TrenaProject", "extdata", "genomeAnnotation", "geneHancer.v4.7.allGenes.RData"))) # 61
dim(tbl.enhancers)

load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))
indices <- match(rownames(mtx), tbl.geneInfo$ensg)
length(indices)
mtx.geneSymbols <- tbl.geneInfo$geneSymbol[indices]
# 
geneSymbols.with.enhancers <- intersect(mtx.geneSymbols, unique(tbl.enhancers$geneSymbol))
#geneSymbols.with.enhancers <- intersect(rownames(mtx), unique(tbl.enhancers$geneSymbol))

printf("geneSymbols.with.enhancers: %d", length(geneSymbols.with.enhancers)) # 14769

ensg.with.enhancers <- subset(tbl.geneInfo, geneSymbol %in% geneSymbols.with.enhancers)$ensg
#length(unique(ensg.with.enhancers))  # 14776
#ensg.without.enhancers <- setdiff(rownames(mtx), ensg.with.enhancers)
#length(ensg.without.enhancers)       # 11560

#goi.test <- c(ensg.with.enhancers[1], ensg.without.enhancers[1])
#goi <- goi.test
goi <- intersect(rownames(mtx), ensg.with.enhancers)
printf("established %d goi", length(goi))
configurationFileRead <- TRUE
tfPrefilterCorrelation=0.1
correlationThreshold=0.1


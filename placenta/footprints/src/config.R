library(TrenaProjectPlacenta)
library(MotifDb)

stopifnot(packageVersion("TrenaProject") >= "0.99.36")
stopifnot(packageVersion("TrenaProjectPlacenta") >= "0.99.17")

OUTPUTDIR <- "2019mar29"
WORKERS <- 20

SOLVERS <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "sqrtLasso")

fp.logDir <- file.path(OUTPUTDIR, "logs.fp")
tfMapping.logDir <- file.path(OUTPUTDIR, "logs.tfMapping")
model.logDir <- file.path(OUTPUTDIR, "logs.model")

trenaProject <- TrenaProjectPlacenta()
available.footprint.databases <- getFootprintDatabaseNames(trenaProject) # "
desired.footprint.databases <- "placenta2"
stopifnot(all(desired.footprint.databases %in% available.footprint.databases))

getExpressionMatrixNames(trenaProject)
matrix.name <- "CombatAdjustedData_MWHONLY_312019_geneSymbols"
stopifnot(matrix.name %in% getExpressionMatrixNames(trenaProject))
mtx <- getExpressionMatrix(trenaProject, matrix.name)

deleters <- grep("^ENSG", rownames(mtx))
if(length(deleters) > 0){
    printf("deleting %d ENSG rows from mtx", length(deleters))
    mtx <- mtx[-deleters,]
    }

print(load(system.file(package="TrenaProject", "extdata", "genomeAnnotation", "geneHancer.v4.7.allGenes.RData"))) # 61
tbl.geneHancer <- tbl.enhancers
dim(tbl.geneHancer)

load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))
indices <- match(rownames(mtx), tbl.geneInfo$ensg)
length(indices)
geneSymbols.in.mtx <- rownames(mtx)
   # eliminate chrX, chrY, chrM genes, for which we have not yet established footprints
geneSymbols.in.mtx <- subset(tbl.geneInfo, geneSymbol %in% rownames(mtx) & chrom %in% sprintf("chr%d", 1:22))$geneSymbol

# 
geneSymbols.with.enhancers <- intersect(geneSymbols.in.mtx, unique(tbl.geneHancer$geneSymbol))
printf("genes with geneHancer annotation: %d/%d", length(geneSymbols.with.enhancers), length(geneSymbols.in.mtx))
head(setdiff(geneSymbols.in.mtx, unique(tbl.geneHancer$geneSymbol)))
goi <- geneSymbols.in.mtx
printf("established %d goi", length(goi))
configurationFileRead <- TRUE
tfPrefilterCorrelation=0.1
correlationThreshold=0.1
tf.pool <- (intersect(trenaSGM::allKnownTFs(identifierType="geneSymbol"), mcols(MotifDb)$geneSymbol))
tf.pool <- intersect(tf.pool, rownames(mtx))
printf("using %d tfs, each with a MotifDb matrix and expression in mtx", length(tf.pool))
use.geneHancer <- TRUE


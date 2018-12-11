library(TrenaProjectGBM)
library(TrenaProject)

stopifnot(packageVersion("TrenaProject") >= "0.99.24")
stopifnot(packageVersion("TrenaProjectGBM") >= "0.99.03")

OUTPUTDIR <- "demo"
SOLVERS <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "sqrtLasso")

fp.logDir <- "logs.fp"
tfMapping.logDir <- "logs.tfMapping"
model.logDir <- "logs.model"

trenaProject <- TrenaProjectGBM()
getExpressionMatrixNames(trenaProject)
matrix.name <- "Scaled_Winsorized_MayoRNAseq_TCX-ENSG"
stopifnot(matrix.name %in% getExpressionMatrixNames(trenaProject))
mtx <- getExpressionMatrix(trenaProject, matrix.name)

stopifnot(grepl("ENSG", rownames(mtx)[1]))

print(load(system.file(package="TrenaProject", "extdata", "geneInfoTable.RData")))

failures <- which(nchar(tbl.geneInfo$hgnc_symbol) == 0)
length(failures)

dim(mtx)
dim(tbl.geneInfo)
all(rownames(mtx) %in% tbl.geneInfo$ensg)
no.transcript.ensgs <- setdiff(rownames(mtx), unique(tbl.geneInfo$ensg))  # [1] 61

print(load(system.file(package="TrenaProject", "extdata", "epigenome", "geneHancer.v4.7.allGenes.RData"))) # 61
dim(tbl.enhancers)
geneSymbols.with.enhancers <- intersect(tbl.geneInfo$geneSymbol, unique(tbl.enhancers$geneSymbol))
printf("geneSymbols.with.enhancers: %d", length(geneSymbols.with.enhancers)) # 15294

ensg.with.enhancers <- subset(tbl.geneInfo, geneSymbol %in% geneSymbols.with.enhancers)$ensg
length(unique(ensg.with.enhancers))  # 15295, no dups
ensg.without.enhancers <- setdiff(rownames(mtx), ensg.with.enhancers)
length(ensg.without.enhancers)       # 1708

length(setdiff(rownames(mtx), ensg.with.enhancers))   # 1708/17003 = 15295 mtx engs have no enhancers

goi.test <- c(ensg.with.enhancers[1], ensg.without.enhancers[1])
goi <- head(rownames(mtx), n=5)
#goi <- "ENSG00000227232"
printf("established %d goi", length(goi))
configurationFileRead <- TRUE

library(TrenaProjectGBM)
library(TrenaProject)

OUTPUTDIR <- "demo"
SOLVERS <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "sqrtLasso")

fp.logDir <- "logs.fp"
tfMapping.logDir <- "logs.tfMapping"
model.logDir <- "logs.model"

trenaProject <- TrenaProjectGBM()
getExpressionMatrixNames(trenaProject)
mtx <- getExpressionMatrix(trenaProject, "tbl.gbm.endogenous")

load(system.file(package="TrenaProject", "extdata", "geneInfoTable.RData"))

failures <- which(nchar(tbl.geneInfo$hgnc_symbol) == 0)
length(failures)

dim(mtx)
dim(tbl.geneInfo)
all(rownames(mtx) %in% tbl.geneInfo$ensg)
no.transcript.ensgs <- setdiff(rownames(mtx), unique(tbl.geneInfo$ensg))  # [1] 61

print(load(system.file(package="TrenaProjectGBM", "extdata", "epigenome", "geneHancer.v4.7.allGenes.RData"))) # 61
dim(tbl.enhancers)
geneSymbols.with.enhancers <- intersect(tbl.geneInfo$geneSymbol, unique(tbl.enhancers$geneSymbol))
printf("geneSymbols.with.enhancers: %d", length(geneSymbols.with.enhancers)) # 15294

ensg.with.enhancers <- subset(tbl.geneInfo, geneSymbol %in% geneSymbols.with.enhancers)$ensg
length(unique(ensg.with.enhancers))  # 15295, no dups
ensg.without.enhancers <- setdiff(rownames(mtx), ensg.with.enhancers)
length(ensg.without.enhancers)       # 1708

length(setdiff(rownames(mtx), ensg.with.enhancers))   # 1708/17003 = 15295 mtx engs have no enhancers

goi.test <- c(ensg.with.enhancers[1], ensg.without.enhancers[1])
goi <- head(rownames(mtx), n=2)
printf("established %d goi", length(goi))
configurationFileRead <- TRUE

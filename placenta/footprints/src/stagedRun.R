# fastRun.R: experiment with ensembl genes
#------------------------------------------------------------------------------
# we need a configuration file (of R commands), specified on the command line
# informing us of how to run this whole genome parallelized script
#------------------------------------------------------------------------------
if(interactive()){
   startGeneIndex <- 1
   endGeneIndex <- 10
   }
if(!interactive()){
   args <- commandArgs(trailingOnly=TRUE)
   stopifnot(length(args) == 2)
   startGeneIndex <- as.integer(args[1])
   endGeneIndex <- as.integer(args[2])
   }
printf("running with genes %d - %d", startGeneIndex, endGeneIndex)
#-----------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(BiocParallel)
library(futile.logger)
library(RPostgreSQL)
library(org.Hs.eg.db)
library(BatchJobs)
#------------------------------------------------------------------------------
stopifnot(packageVersion("trena") >= "1.5.13")
stopifnot(packageVersion("trenaSGM") >= "0.99.76")
#------------------------------------------------------------------------------
configurationFile <- "config.R"
stopifnot(file.exists(configurationFile))

if(!exists("configurationFileRead"))
  source(configurationFile)

stopifnot(exists("trenaProject"))
stopifnot(exists("mtx"))
stopifnot(exists("goi"))


if(!file.exists(OUTPUTDIR)) dir.create(OUTPUTDIR)
if(!file.exists(fp.logDir)) dir.create(fp.logDir)
if(!file.exists(tfMapping.logDir)) dir.create(tfMapping.logDir)
if(!file.exists(model.logDir)) dir.create(model.logDir)

#----------------------------------------------------------------------------------------------------
basic.build.spec <- list(title="footprint-based-tf-model-builder-for-placental-genes",
                         type="footprint.database",
                         stageDirectory=OUTPUTDIR,
                         genomeName="hg38",
                         matrix=mtx,
                         db.host=getFootprintDatabaseHost(trenaProject),
                         databases=desired.footprint.databases,
                         annotationDbFile=dbfile(org.Hs.eg.db),
                         motifDiscovery="builtinFimo",
                         tfPool=tf.pool,
                         tfMapping=c("MotifDB", "TFClass"),
                         tfPrefilterCorrelation=tfPrefilterCorrelation, 
                         correlationThreshold=correlationThreshold,
                         orderModelByColumn="rfScore",
                         solverNames=SOLVERS,
                         solvers=SOLVERS)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_determineRegulatoryRegions()

} # runTests
#----------------------------------------------------------------------------------------------------
# cory's method:
# if there is no ENSG in tbl.geneInfo, return an empty data.frame
# if there is no geneSymbol for a given ENSG id, then no enhancers are avaialable: use tss +/- 5k
# if there is a geneSymbol, but no enhancers entry: use tss +/- 5k
# if the ENSG has a geneSymbol, and the geneSymbol has associated enhancers: use those enhancers 
#
# --  alison's strategy (email 26 mar 2019)
#
#  I propose that we build this first model using +/- 2500 BP from the TSS & any distal enhancer regions
#  outside of this region identified within genehancer. I think that this strikes a good balance of being as
#  inclusive of a broad range of the promoter as possible, and also will capture distal hits in a way that is
#  fair. This was the most liberal of the 3 strategies I came up with. We can tone this down if the models end
#  up looking strange.  Let me know if this makes sense and what you think of this strategy. Thanks!

determineRegulatoryRegions <- function(gene)
{
   shoulder <- 2500
   tbl.concise <- tbl.geneInfo[grep(gene, tbl.geneInfo$geneSymbol), c("chrom", "tss")]
      # no need to figure strand since we go 2500bp in both directions
   tbl.regions <- data.frame(chrom=tbl.concise$chrom,
                             start=tbl.concise$tss - shoulder,
                             end=tbl.concise$tss + shoulder,
                             stringsAsFactors=FALSE)
   
   setTargetGene(trenaProject, gene)
   tbl.enhancers <- data.frame()
   if(use.geneHancer)
      tbl.enhancers <- getEnhancers(trenaProject)[, c("chrom", "start", "end")]

   tbl.regions <- rbind(tbl.regions, tbl.enhancers)
   
   new.order <- order(tbl.regions$start)
   tbl.regions <- tbl.regions[new.order,]
   rownames(tbl.regions) <- NULL
   tbl.reduced <- as.data.frame(union(GRanges(tbl.regions), GRanges(tbl.regions)))[, c("seqnames", "start", "end")]
   colnames(tbl.reduced) <- c("chrom", "start", "end")
   tbl.reduced$chrom <- as.character(tbl.reduced$chrom)

   return(tbl.reduced)

} # determineRegulatoryRegions
#----------------------------------------------------------------------------------------------------
test_determineRegulatoryRegions <- function()
{
   printf("--- test_determineRegulatoryRegions")

   tbl.tspan6 <- determineRegulatoryRegions("TSPAN6")
   checkEquals(dim(tbl.tspan6), c(3, 3))
   tbl.tut7 <- determineRegulatoryRegions("TUT7")
   checkEquals(dim(tbl.tut7), c(1,3))

} # test_determineRegulatoryRegions
#----------------------------------------------------------------------------------------------------
runStagedSGM.footprints <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol") # , "regionsMode")
   missing.fields <- setdiff(required.fields, names(short.spec))

   if(length(missing.fields) > 0){
      msg <- sprintf("runStagedSGM.footprints finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   printf("-- runSGM(%s:%s) footprints", short.spec$targetGene, short.spec$geneSymbol)

   targetGene <- short.spec$targetGene
   geneSymbol <- short.spec$geneSymbol

   stopifnot(nchar(targetGene) >= 2)
   stopifnot(nchar(geneSymbol) >= 2)

   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss

   tbl.regions <- determineRegulatoryRegions(targetGene) # , short.spec$regionsMode)
   

   spec <- basic.build.spec
   spec$regions <- tbl.regions
   spec$targetGene <- targetGene
   spec$tss=tss
   spec$geneSymbol=geneSymbol
   printf("--- path: %s:", spec$stageDir)

   fpBuilder <- FastFootprintDatabaseModelBuilder(spec$genomeName, targetGene, spec,
                                                  quiet=FALSE,
                                                  stagedExecutionDirectory=spec$stageDir)
   fp.filename <- staged.fast.build(fpBuilder, stage="find.footprints")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.footprints
#----------------------------------------------------------------------------------------------------
runStagedSGM.associateTFs <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol") #, "regionsMode", "correlationThreshold", "solvers", "dbs")
   missing.fields <- setdiff(required.fields, names(short.spec))
   if(length(missing.fields) > 0){
      msg <- sprintf("runStagedSGM.footprints finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   printf("-- runSGM, assoicateTFs(%s)", short.spec$targetGene)

   genomeName <- "hg38"
   targetGene <- short.spec$targetGene
   geneSymbol <- short.spec$geneSymbol
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss

   spec <- basic.build.spec
   spec$targetGene <- targetGene
   spec$tss=tss
   spec$regions <- data.frame()
   spec$geneSymbol=geneSymbol

   printf("--- path: %s:", spec$stageDir)

   fpBuilder <- FastFootprintDatabaseModelBuilder(genomeName, targetGene, spec, quiet=FALSE,
                                                  stagedExecutionDirectory=OUTPUTDIR)
   fp.filename <- staged.fast.build(fpBuilder, stage="associateTFs")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.associateTFs
#----------------------------------------------------------------------------------------------------
runStagedSGM.buildModels <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol", "regionsMode") # , "correlationThreshold", "solvers", "dbs")
   missing.fields <- setdiff(required.fields, names(short.spec))
   if(length(missing.fields) > 0){
      msg <- sprintf("runStagedSGM.buildModels finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   printf("-- runSGM(%s) buildModels", short.spec$targetGene)

   genomeName <- "hg38"
   targetGene <- short.spec$targetGene
   geneSymbol <- short.spec$geneSymbol
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss

   #tbl.regions <- switch(short.spec$regionsMode,
   #                      "enhancers" = {enhancer.list[[targetGene]][, c("chrom", "start", "end")]},
   #                      "tiny" = {data.frame(chrom=chromosome, start=tss-1000, end=tss+1000, stringsAsFactors=FALSE)
   #                      })

   spec <- basic.build.spec
   spec$targetGene <- targetGene
   spec$tss=tss
   spec$regions <- data.frame()
   spec$geneSymbol=geneSymbol
   printf("--- path: %s:", spec$stageDir)

   fpBuilder <- FastFootprintDatabaseModelBuilder(genomeName, targetGene, spec, quiet=FALSE,
                                              stagedExecutionDirectory=OUTPUTDIR)
   fp.filename <- staged.fast.build(fpBuilder, stage="build.models")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.buildModels
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.footprints <- function(genes)
{
   printf("--- kicking off do.runStagedSGM.footprints, gene count %d", length(genes))
    
   short.specs <- lapply(genes,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=gene,
                                 regionsMode="enhancers"))
   names(short.specs) <- as.character(genes)

   if(interactive()) runStagedSGM.footprints(short.specs[[1]])   # WASH7P
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=fp.logDir, threshold="INFO", workers=WORKERS)
   results.fp <<- bptry({bplapply(short.specs, runStagedSGM.footprints, BPPARAM=bp.params)})

} # do.runStagedSGM.footprints
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.associateTFs <- function(genes)
{
   printf("--- kicking off do.runStagedSGM.associateTFs, gene count %d", length(genes))

   short.specs <- lapply(genes,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=gene,
                                 regionsMode="enhancers"))
   names(short.specs) <- as.character(genes)

   if(interactive()) runStagedSGM.associateTFs(short.specs[[1]])
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=tfMapping.logDir, threshold="INFO", workers=WORKERS)
   results.tfMap <<- bptry({bplapply(short.specs, runStagedSGM.associateTFs, BPPARAM=bp.params)})

} # do.runStagedSGM.associateTFs
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.buildModels <- function(genes)
{
   printf("--- kicking off do.runStagedSGM.buildModels, gene count %d", length(genes))

   short.specs <- lapply(genes,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=gene,
                                 regionsMode="enhancers"))

   names(short.specs) <- as.character(genes)

   if(interactive()) runStagedSGM.buildModels(short.specs[[1]])

   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=model.logDir, threshold="INFO", workers=WORKERS)
   results.buildModels <<- bptry({bplapply(short.specs, runStagedSGM.buildModels, BPPARAM=bp.params)})

} # do.runStagedSGM.buildModels
#----------------------------------------------------------------------------------------------------
genes <- goi[startGeneIndex:endGeneIndex]
genes <- intersect(genes, tbl.geneInfo$geneSymbol)
printf("total genes this run: %d", length(genes))

printf("starting run of %d (%d-%d) genes, writing staged results to %s", length(genes),
        startGeneIndex, endGeneIndex, OUTPUTDIR)

do.runStagedSGM.footprints(genes)
do.runStagedSGM.associateTFs(genes)
do.runStagedSGM.buildModels(genes)

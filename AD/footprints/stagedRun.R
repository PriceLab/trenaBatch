# fastRun.R: experiment with ensembl genes
#------------------------------------------------------------------------------
# we need a configuration file (of R commands), specified on the command line
# informing us of how to run this whole genome parallelized script
#------------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(motifStack)
library(tibble)
library(BiocParallel)
library(futile.logger)
library(RPostgreSQL)
library(org.Hs.eg.db)
library(BatchJobs)
#------------------------------------------------------------------------------
stopifnot(packageVersion("TrenaProject") >= "0.99.22")
stopifnot(packageVersion("TrenaProjectIGAP") >= "0.99.14")
stopifnot(packageVersion("trena") >= "1.5.1")
stopifnot(packageVersion("trenaSGM") >= "0.99.57")
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
basic.build.spec <- list(title="footprint-based-tf-model-builder-for-ADgenes",
                         type="footprint.database",
                         stageDirectory=OUTPUTDIR,
                         genomeName="hg38",
                         matrix=mtx,
                         db.host=getFootprintDatabaseHost(trenaProject),
                         databases=getFootprintDatabaseNames(trenaProject)[1],
                         annotationDbFile=dbfile(org.Hs.eg.db),
                         motifDiscovery="builtinFimo",
                         tfPool=allKnownTFs(identifierType="ensemblGeneID"),
                         tfMapping=c("MotifDB", "TFClass"),
                         tfPrefilterCorrelation=0.2,
                         correlationThreshold=0.2,
                         orderModelByColumn="rfScore",
                         solverNames=SOLVERS,
                         solvers=SOLVERS)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_determineRegulatoryRegions()

} # runTests
#----------------------------------------------------------------------------------------------------
# if(!exists("targetGenes")){   # a named list, geneSymbols as names, ensg as content
#    indices <- match(goi, rownames(tbl.geneInfo))
#    deleters <- which(is.na(indices))
#    if(length(deleters) > 0){
#       goi <- goi[-deleters]
#       indices <- indices[-deleters]
#       }
#    targetGenes <- rownames(tbl.geneInfo[indices,])
#    names(targetGenes) <- goi
#    }
#----------------------------------------------------------------------------------------------------
# cory's method:
# if there is no ENSG in tbl.geneInfo, return an empty data.frame
# if there is no geneSymbol for a given ENSG id, then no enhancers are avaialable: use tss +/- 5k
# if there is a geneSymbol, but no enhancers entry: use tss +/- 5k
# if the ENSG has a geneSymbol, and the geneSymbol has associated enhancers: use those enhancers AND +/- 2k
determineRegulatoryRegions <- function(gene)
{
   BIG.SHOULDER <- 5000     # used when no enhancers are available
   SMALL.SHOULDER <- 1000   # used to complement enhancers

   stopifnot(grepl("^ENSG", gene))

     # need geneSymbol to look up enhancers
   tbl.regulatoryRegions.failed <- data.frame()

   index <- match(gene, tbl.geneInfo$ensg)
   if(all(is.na(index)))
      return(tbl.regulatoryRegions.failed)

   index <- index[1]   # just in case
   tbl.thisGene <- tbl.geneInfo[index,]
   tbl.out <- with(tbl.thisGene, data.frame(chrom=chrom,
                                            start=tss-BIG.SHOULDER,
                                            end=tss+BIG.SHOULDER,
                                            stringsAsFactors=FALSE))

   geneSymbol <- tbl.thisGene$geneSymbol
   tbl.enhancers <- getEnhancers(trenaProject, geneSymbol)

   if(nrow(tbl.enhancers) == 0){
      return(tbl.out)
      }

   tbl.out <- with(tbl.thisGene, data.frame(chrom=chrom,
                                            start=tss-SMALL.SHOULDER,
                                            end=tss+SMALL.SHOULDER,
                                            stringsAsFactors=FALSE))

   tbl.enhancers <- tbl.enhancers[,c("chrom", "start", "end")]

   tbl.out <- rbind(tbl.out, tbl.enhancers)
   new.order <- order(tbl.out$start)
   tbl.out <- tbl.out[new.order,]
   rownames(tbl.out) <- NULL

   return(unique(tbl.out))

} # determineRegulatoryRegions
#----------------------------------------------------------------------------------------------------
test_determineRegulatoryRegions <- function()
{
   printf("--- test_determineRegulatoryRegions")

      # no info on this bogus gene:
   checkEquals(nrow(determineRegulatoryRegions("ENSGhocusPocusGene")), 0)

      # only tss, no enhancers
   tbl.1 <- determineRegulatoryRegions("ENSG00000279457")
   checkEquals(nrow(tbl.1), 1)
   checkEquals(with(tbl.1, end-start), 10000)

     # now an ENSG gene for which we have tss AND enhancer info
   tbl.regions <- determineRegulatoryRegions("ENSG00000000003")
   checkEquals(dim(tbl.regions), c(4,3))
   checkEquals(colnames(tbl.regions), c("chrom", "start", "end"))
   checkEquals(with(tbl.regions, end-start), c(1056, 2174, 2000, 8909))
   checkEquals(with(tbl.regions, sum(end-start)), 14139)

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

   printf("-- runSGM(%s) footprints", short.spec$targetGene)

   targetGene <- short.spec$targetGene
   geneSymbol <- short.spec$geneSymbol
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
do.runStagedSGM.footprints <- function()
{
   short.specs <- lapply(goi,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=subset(tbl.geneInfo, ensg==gene)$geneSymbol,
                                 regionsMode="enhancers"))
   names(short.specs) <- as.character(goi)

   if(interactive()) runStagedSGM.footprints(short.specs[[1]])   # WASH7P
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=fp.logDir, level="INFO")
   results.fp <<- bptry({bplapply(short.specs, runStagedSGM.footprints, BPPARAM=bp.params)})

} # do.runStagedSGM.footprints
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.associateTFs <- function()
{
   short.specs <- lapply(goi,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=subset(tbl.geneInfo, ensg==gene)$geneSymbol,
                                 regionsMode="enhancers"))
   names(short.specs) <- as.character(goi)

   if(interactive()) runStagedSGM.associateTFs(short.specs[[1]])
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=tfMapping.logDir, level="INFO")
   results.tfMap <<- bptry({bplapply(short.specs, runStagedSGM.associateTFs, BPPARAM=bp.params)})

} # do.runStagedSGM.associateTFs
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.buildModels <- function()
{
   #ensembl.tfPool <- allKnownTFs(identifierType="ensemblGeneID")

   short.specs <- lapply(goi,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=subset(tbl.geneInfo, ensg==gene)$geneSymbol,
                                 regionsMode="enhancers"))

   names(short.specs) <- as.character(goi)

   if(interactive()) runStagedSGM.buildModels(short.specs[[1]])

   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=model.logDir, level="INFO")
   results.buildModels <<- bptry({bplapply(short.specs, runStagedSGM.buildModels, BPPARAM=bp.params)})

} # do.runStagedSGM.buildModels
#----------------------------------------------------------------------------------------------------
if(!interactive()){
   printf("starting run of %d goi, writing staged results to %s", length(goi), OUTPUTDIR)
   do.runStagedSGM.footprints()
   do.runStagedSGM.associateTFs()
   do.runStagedSGM.buildModels()
   }

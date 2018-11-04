library(BiocParallel)
library(trenaSGM)
library(TrenaProjectPlacenta)
library(futile.logger)
#------------------------------------------------------------------------------------------------------------------------
project <- TrenaProjectPlacenta()
mtx <- getExpressionMatrix(project, "FilteredCountData8282018-vsn-geneSymbols")
deleters <- which(is.na(rownames(mtx)))
if(length(deleters) > 0)
    mtx <- mtx[-deleters,]
dim(mtx)
goi <- getSupportedGenes(project)
#------------------------------------------------------------------------------------------------------------------------
runGeneModel <- function(project, gene)
{
  tryCatch({
  printf("--- runGeneModel %s", gene)
  filename.out <- sprintf("out/%s.out.RData", gene)
  printf("filename.out:")
  print(filename.out)
  if(file.exists(filename.out)){
      printf("%s already exists, skipping", filename.out)
      return()
      }

  setTargetGene(project, gene)
  tbl.transcripts <- getTranscriptsTable(project)
  if(nrow(tbl.transcripts) == 0){
      printf("no tbl.transcripts for %s, skipping", gene)
      return()
      }
      
  tbl.geneData <- subset(tbl.transcripts, moleculetype=="gene")[1,]

  chrom <- tbl.geneData$chr
  start <- tbl.geneData$start
  end  <- tbl.geneData$endpos
  strand <- tbl.geneData$strand
  tss <- start
  if(strand == "-")
     tss <- end

  tbl.enhancers <- getEnhancers(project)
  if(nrow(tbl.enhancers) > 0){ 
     tbl.enhancers <- tbl.enhancers[order(tbl.enhancers$combinedScore, decreasing=TRUE),]
     regions.to.use <- 3
     if(nrow(tbl.enhancers) < 3)
         regions.to.use <- nrow(tbl.enhancers)
     tbl.regions <- tbl.enhancers[1:regions.to.use,]
    }
  else{
      tbl.regions <- data.frame(chrom=chrom, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)
      }
  
  recipe <- list(title=sprintf("motifMatching.%s", getTargetGene(project)),
                 type="regions.motifMatching",
                  geneSymbol=gene,
                  tss=tss,
                  regions=tbl.regions,
                  matrix=mtx,
                  pfms=query(MotifDb, "sapiens", "jaspar2018"),
                  matchThreshold=80,
                  motifDiscovery="matchPWM",
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.1,
                  orderModelByColumn="rfScore",
                  annotationDbFile=dbfile(org.Hs.eg.db),
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                  quiet=FALSE)

    builder <- RegionsMotifMatchingModelBuilder("hg38", gene, recipe, quiet=FALSE)
    times <- as.list(system.time(x <- build(builder)))
    printf("%s time %f", gene, times$elapsed)
    save(x, file=filename.out)
  }, error=function(e) {
      printf("failure with %s: %s", gene, e)
      return(list(model=data.frame(), regulatoryRegions=data.frame()))
      })


} # runGeneModel
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()){
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir="out/logs", level="INFO", workers=20)
   #bp.parmas <- MulticoreParam(workers=10)
   x <- list.files("out")
   gad <- sub(".out.RData", "", x)
   goi <- setdiff(rownames(mtx), gad)
   printf("%d genes remaining to model", length(goi))
   lapply(goi, function(gene) runGeneModel(project, gene))
   #bplapply(rownames(mtx), function(gene) runGeneModel(project, gene), BPPARAM=bp.params)

   }



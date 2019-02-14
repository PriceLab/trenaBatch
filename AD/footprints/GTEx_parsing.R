# large matrix and other files are located on whovian:/local/Cory/gtex

# load the matrix
mtx <- readRDS("GeneReads.rds")

# print(load("matrix.colnames.RData"))
# length(matrix.colnames)
# head(matrix.colnames)
tbl.md <- read.table("GTEx_v7_Annotations_SampleAttributesDS.txt", sep = '\t', header=TRUE,
                     quote="", as.is=TRUE, row.names=1)

dim(tbl.md) # 15598 63

matrix.colnames <- colnames(mtx)
length(intersect(matrix.colnames, rownames(tbl.md)))

new.colnames <- gsub(".", "-", matrix.colnames, fixed=TRUE)
length(intersect(new.colnames, rownames(tbl.md)))

colnames(mtx) <- new.colnames

tbl.tissues <- as.data.frame(table(tbl.md$SMTSD), stringsAsFactors=FALSE)
colnames(tbl.tissues) <- c("tissue", "sampleCount")
descending.order <- order(tbl.tissues$sampleCount, decreasing=TRUE)
tbl.tissues <- tbl.tissues[descending.order,]

#  identify columns (samples) for a tissue with only a few samples
tail(tbl.tissues)

# get the Whole_blood samples (total of 2412)
blood.sampleIDs <- rownames(subset(tbl.md, SMTSD=="Whole_Blood"))
length(blood.sampleIDs)
length(intersect(blood.sampleIDs, colnames(mtx)))

length(blood.sampleIDs) # 2412

blood.intersect <- which(blood.sampleIDs %in% colnames(mtx))
length(blood.intersect) # 407

mtx.tissue <- mtx[, blood.intersect]

# get rid of the decimal points in the ensembl ids
rownames(mtx.tissue) <- sapply(strsplit(rownames(mtx.tissue),"\\."), `[`, 1)
#----------------------------------------------------------
# loop through all the tissues and check the sample numbers
described.count <- 0
measured.count <- 0

for(tissueName in tbl.tissues$tissue){
   samples.described <- rownames(subset(tbl.md, SMTSD==tissueName))
   samples.measured <- length(intersect(samples.described, colnames(mtx)))
   described.count <- described.count + length(samples.described)
   measured.count <- measured.count + samples.measured
   printf("tissueName: %s: %d/%d", tissueName, samples.measured, length(samples.described))
   }

printf("described: %d   measured:  %d", described.count, measured.count)
#----------------------------------------------------------


# extract from the big all-tissue matrix
#mtx.blood <- mtx[, blood.sampleIDs]

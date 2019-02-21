target.dir <- "2019feb15"
load(file.path(target.dir, "tbl.models.all.RData"))
dim(tbl.models)

tbl.10 <- subset(tbl.models, rank <= 10)
dim(tbl.10)
tbl.tf10.dist <- as.data.frame(table(tbl.10$tf.symbol))
colnames(tbl.tf10.dist) <- c("tf", "count")
tbl.tf10.dist <- tbl.tf10.dist[order(tbl.tf10.dist$count, decreasing=TRUE),]
head(tbl.tf10.dist)

tbl.05 <- subset(tbl.models, rank <= 05)
dim(tbl.05)
tbl.tf05.dist <- as.data.frame(table(tbl.05$tf.symbol))
colnames(tbl.tf05.dist) <- c("tf", "count")
tbl.tf05.dist <- tbl.tf05.dist[order(tbl.tf05.dist$count, decreasing=TRUE),]
head(tbl.tf05.dist)

fivenum(tbl.tf10.dist$count)
fivenum(tbl.tf05.dist$count)

fivenum(abs(tbl.10$spearmanCoeff))   # [1] 0.0457974 0.4338310 0.5399962 0.6431475 0.9586307
fivenum(abs(tbl.05$spearmanCoeff))   # [1] 0.0457974 0.4667007 0.5709940 0.6707557 0.9586307

hist(tbl.05$spearmanCoeff)
hist(tbl.10$spearmanCoeff)

hist(tbl.10$rfScore)
hist(tbl.05$rfScore)


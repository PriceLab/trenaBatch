max <- 13784
starts <- seq(1, 13800, 20)
ends <- starts + 19
ends[length(ends)] <- max
bashScript <- file("runByChunks.sh")
lines <- sprintf("Rscript stagedRun.R %d %d", starts, ends)
writeLines(lines, bashScript)
close(bashScript)

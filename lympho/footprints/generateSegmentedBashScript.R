max <- 26172
starts <- seq(1, 27200, 100)
ends <- starts + 99
ends[length(ends)] <- max
bashScript <- file("runBy100s.sh")
lines <- sprintf("Rscript ../../stagedRun.R %d %d", starts, ends)
writeLines(lines, bashScript)
close(bashScript)

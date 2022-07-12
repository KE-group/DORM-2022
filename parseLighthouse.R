rm(list = ls())
gc()
library(stringi)
# library(data.table)

setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/Performance Comparison/Lighthouse reports/")

# grep JSON from lighthouse reports using ripgrep.
system("rg '<script>window.__LIGHTHOUSE_JSON__ =' >| ripgrep_output.txt", intern = F, wait = F)
Data <- readLines(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/Performance Comparison/Lighthouse reports/ripgrep_output.txt"
)

to_keep <-
  c(
    "firstContentfulPaint",
    "firstMeaningfulPaint",
    "largestContentfulPaint",
    "interactive",
    "speedIndex",
    "totalBlockingTime",
    "cumulativeLayoutShift"
  )
summaryDF <-
  data.frame(matrix(NA, nrow = length(Data), ncol = (length(to_keep) + 2)))
colnames(summaryDF) <- c("filename", "performanceScore", to_keep)

for (i in seq(1, length(Data))) {
  LINE <- Data[i]
  
  var <- stringi::stri_split_fixed(str = LINE, pattern = ".html: ")
  filename <- gsub("/", "_", var[[1]][1])
  JSON <- var[[1]][2]
  rm(var)
  gc()
  
  # Cleaning up JSON
  JSON <- stringi::stri_replace_first_fixed(str = JSON,
                                            pattern = "<script>window.__LIGHTHOUSE_JSON__ = ",
                                            replacement = "")
  JSON <- stringi::stri_replace_first_fixed(str = JSON,
                                            pattern = ";</script>",
                                            replacement = "")
  
  # Parsing JSON
  raw_json <- jsonlite::fromJSON(txt = JSON)
  
  tmp <-
    t(raw_json[["audits"]][["metrics"]][["details"]][["items"]])
  metrics <- data.frame(name = rownames(tmp), value = tmp[, 1])
  rownames(metrics) <- NULL
  
  metrics[metrics$name %in% to_keep, "value"]
  summaryDF[i, ] <- c(filename,
                      raw_json[["categories"]][["performance"]][["score"]],
                      metrics[metrics$name %in% to_keep, "value"])
  
}

summaryDF <- summaryDF[order(summaryDF$filename, decreasing = F), ]
write.table(
  x = summaryDF,
  file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/Performance Comparison/LighthouseStats.csv",
  quote = F,
  sep = ",",
  row.names = F,
  col.names = T
)

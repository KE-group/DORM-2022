#------> data.frame vs data.table ------
# Count the total number of samples from each tissue
test.name <- "grep"

Data <-
  readLines(
    "/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/websites/eleniuslabtools.utu.fi/shiny-apps/DORM/Mutations/data/tissue/all.csv"
  )

bmark_out_tsv <- paste0("table/bmark", "_", test.name, ".tsv")
file.create(bmark_out_tsv, showWarnings = F)

searchExpression <- "EGFR|ERBB[2-4]|[HKN]RAS\\>"

GNU_grep_command <-
  paste0(
    "/usr/local/opt/grep/bin/gegrep -i '",
    searchExpression,
    "' '/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/websites/eleniuslabtools.utu.fi/shiny-apps/DORM/Mutations/data/tissue/all.csv' > /dev/null"
  )

BSD_grep_command <-
  paste0(
    "/usr/bin/egrep -i '",
    searchExpression,
    "' '/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/websites/eleniuslabtools.utu.fi/shiny-apps/DORM/Mutations/data/tissue/all.csv' > /dev/null"
  )


bmark <- microbenchmark(
  "base_r" = {
    write(x = Data[grep(pattern = searchExpression,
                        x = Data,
                        ignore.case = T)],
          file = "/dev/null",
          ncolumns = 1)
  },
  "stringi" = {
    write(x = Data[stringi::stri_detect_regex(str = Data,
                                              pattern = searchExpression,
                                              case_insensitive = T)],
          file = "/dev/null",
          ncolumns = 1)
  },
  "BSD_grep" = {
    system(command = BSD_grep_command,
           intern = F,
           wait = T)
  },
  "GNU_grep" = {
    system(command = GNU_grep_command,
           intern = F,
           wait = T)
  },
  times = 10
)
# saveRDS(bmark, file = paste0("./bmark/bmark_", test.name, "_", i, ".RDS"))
write.table(
  bmark,
  file = bmark_out_tsv,
  sep = "\t",
  row.names = F,
  col.names = F,
  append = T
)

rm(list = ls())
gc()
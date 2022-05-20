IBMdataImport <- read.table(file.path(getwd(), "data-raw/IBM.csv"),
                      quote="\"", comment.char="")
IBMdata <- data.frame(IBMdataImport[2013:5283, ])
names(IBMdata) <- "r_t"
# check if everything looks good:
tibble::tibble(IBMdata)
usethis::use_data(IBMdata)

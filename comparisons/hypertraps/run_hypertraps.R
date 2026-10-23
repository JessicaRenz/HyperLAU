library(Rcpp)
sourceCpp("hypertraps-r.cpp")

data <- read.table(
  "tb_data_10_train_qm30_hypertraps.txt",
  colClasses = c("character", "character")
)

toBinaryMatrix <- function(x) {
  matrix(
    as.numeric(unlist(strsplit(x, ""))),
    ncol = nchar(x[1]),
    byrow = TRUE
  )
}

d_a <- toBinaryMatrix(data[,1])
d_d <- toBinaryMatrix(data[,2])


start_time <- Sys.time()
res=HyperTraPS(obs = d_d, initialstates = d_a)

end_time <- Sys.time()
runtime <- end_time - start_time

outfile <- "hypertraps_output_10_qm30.txt"

cat(
  paste0(
    "HyperTraPS results\n",
    "================\n\n",
    "Runtime: ", runtime, "\n\n"
  ),
  file = outfile
)

# Transitions
cat("\nTRANSITIONS\n-----------\n",
    file = outfile, append = TRUE)

write.table(
  res$edges,
  file = outfile,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  append = TRUE
)


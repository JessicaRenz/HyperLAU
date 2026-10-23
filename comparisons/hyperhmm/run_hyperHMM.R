library(Rcpp)
library(RcppArmadillo)
sourceCpp("hyperhmm-r.cpp")

data <- read.table(
  "tb_data_3_train.txt",
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


same <- apply(d_a == d_d, 1, all)
which(same)

keep <- !same

d_a_red <- d_a[keep, , drop = FALSE]
d_d_red <- d_d[keep, , drop = FALSE]


start_time <- Sys.time()
res = HyperHMM(d_d_red, initialstates=d_a_red,nboot=0)

end_time <- Sys.time()
runtime <- end_time - start_time


outfile <- "hyperhmm_output_3.txt"

# Header + Laufzeit
cat(
  paste0(
    "HyperHMM results\n",
    "================\n\n",
    "Runtime: ", runtime, "\n\n"
  ),
  file = outfile
)

# Transitions
cat("\nTRANSITIONS\n-----------\n",
    file = outfile, append = TRUE)

write.table(
  res$transitions,
  file = outfile,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  append = TRUE
)

# Stats
cat("\n\nSTATS\n-----\n",
    file = outfile,
    append = TRUE)

write.table(
  res$stats,
  file = outfile,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  append = TRUE
)

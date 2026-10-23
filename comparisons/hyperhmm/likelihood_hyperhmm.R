
#provide string length
L = 10
#read in data
data <- read.table(
  "tb_data_10_train.txt",
  stringsAsFactors = FALSE,
  colClasses = c("character","character")
)

colnames(data) <- c("ancestor", "descendant")

#integer to binary
int_to_bits <- function(x, L) {
  
  bits <- as.integer(intToBits(x))[1:L]
  
  rev(bits)
}

#rewrite dataset
data$ancestor_id <- sapply(data$ancestor, state_to_int)
data$descendant_id <- sapply(data$descendant, state_to_int)

is_possible_pair <- function(a, d){
  
  a_bits <- as.integer(strsplit(as.character(a), "")[[1]])
  d_bits <- as.integer(strsplit(as.character(d), "")[[1]])
  
  all(a_bits <= d_bits)
}

keep <- mapply(
  is_possible_pair,
  data$ancestor,
  data$descendant
)

data <- data[keep, ]

#read in hyperHMM output
lines <- readLines("hyperhmm_output_10.txt")

#find transitions block
start <- grep("^TRANSITIONS$", lines)
end <- grep("^STATS$", lines)

trans_lines <- lines[(start + 2):(end - 3)]

transitions <- read.table(
  text = trans_lines,
  header = TRUE,
  sep = "\t"
)

n_states = 2^L

Tmat <- matrix(
  0,
  nrow = n_states,
  ncol = n_states
)

#create transition matrix
for(i in seq_len(nrow(transitions))){
  
  from <- transitions$From[i] + 1
  to <- transitions$To[i] + 1
  
  Tmat[to, from] <- transitions$Probability[i]
}


P_desc <- matrix(
  0,
  nrow = n_states,
  ncol = L + 1
)

P_desc[,1] <- p0

for(t in 2:(L+1)){
  P_desc[,t] <- Tmat %*% P_desc[,t-1]
}

T_powers <- vector("list", L + 1)

T_powers[[1]] <- diag(n_states)

for(k in 1:L){
  T_powers[[k+1]] <- Tmat %*% T_powers[[k]]
}
pair_likelihood <- function(a, d,
                            P_desc,
                            T_powers,
                            L){
  
  lh <- 0
  
  for(t in 0:L){
    
    p_a <- P_desc[a + 1, t + 1]
    
    if(p_a == 0)
      next
    
    for(s in t:L){
      
      delta <- s - t
      
      p_d_given_a <-
        T_powers[[delta + 1]][d + 1, a + 1]
      
      lh <- lh + p_a * p_d_given_a
    }
  }
  
  factor <- (L + 1) * (L + 2)
  
  2 * lh / factor
}

pair_probs <- mapply(
  pair_likelihood,
  a = data$ancestor_id,
  d = data$descendant_id,
  MoreArgs = list(
    P_desc = P_desc,
    T_powers = T_powers,
    L = L
  )
)


sum(log(pair_probs))



### -------- helper functions ------- 

### todo: update functions with docline and parameter input assumptions


## ----- functions to apply to time series data ------ 

### function to split a time series into L blocks of size n = T/L
make_blocks <- function(Tlen, L) {
  stopifnot(Tlen %% L == 0)
  n <- Tlen / L
  split(seq_len(Tlen), rep(seq_len(L), each = n))
}

### fetch ell'th (l) block from the L blocks - corresponds to evaluation data
II_l <- function(Tlen, L, l) { 
  stopifnot(1 <= l & l <= L)
  make_blocks(Tlen, L)[[l]]
}

### fetch l'th training data according to equation 
II_bar_l <- function(Tlen, L, l, p) {
  idx_eval <- II_l(Tlen, L, l)
  
  if (L == 1) { #if no cross-fitting, return whole dataset
    return(seq_len(Tlen))
  }
  
  left  <- max(1, min(idx_eval) - p)
  right <- min(Tlen, max(idx_eval) + p)
  
  setdiff(seq_len(Tlen), left:right)
}

### get pairs (t, t + 1) such that both entries are training points
pairs_from_mask <- function(in_mask) {
  Tlen <- length(in_mask)
  which(in_mask[-Tlen] & in_mask[-1])
}

#### Gets Y as a matrix
get_Y_mat <- function(data) {
  nm <- names(data)
  ycols <- grep("^Y(\\d+)?$", nm)
  if (length(ycols) == 0) stop("No Y columns found. Provide 'Y' or 'Y1..Yd'.")
  as.matrix(data[, ycols, drop = FALSE])
}
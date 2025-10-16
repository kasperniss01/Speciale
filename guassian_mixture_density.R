gaussian_mixture_density <- function(x, weights, means, sds) {
  stopifnot(length(weights) == length(means),
            length(means)   == length(sds))
  
  # Precompute column constants once
  col_const <- weights / (sds * sqrt(2 * pi))
  
  # N x G differences, normalized by sds (column-wise)
  Z <- sweep(outer(x, means, "-"), 2, sds, "/")
  
  # N x G densities (without the column constants)
  dens <- exp(-0.5 * Z^2)
  
  # Weighted sum across components (G) -> length-N vector
  drop(dens %*% col_const)
}


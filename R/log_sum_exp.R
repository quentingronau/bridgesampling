## Log-scale arithmetic helpers.
## These replace the Brobdingnag "brob" class: all quantities stay on the
## log scale until a final exp(), so that very large or very small values
## do not overflow or underflow.

## Elementwise log(exp(a) + exp(b)).
.log_add_exp <- function(a, b) {
  m <- pmax(a, b)
  d <- abs(a - b)
  ## a == b also covers a = b = Inf and a = b = -Inf, where a - b is NaN
  d[which(a == b)] <- 0
  m + log1p(exp(-d))
}

## log(sum(exp(x))) for a vector x.
.log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

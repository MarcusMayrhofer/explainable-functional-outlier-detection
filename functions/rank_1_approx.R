######
# Computation of the best rank-1 Kronecker approximation. This is only stored for the simulations but not used in the paper
# Code by Chat GPT 10.10.2025
######


# Corrected Van-Loan & Pitsianis best rank-1 Kronecker approximation 
# Conventions:
#  - vec(M)  == as.vector(M) (R's column-major stacking)
#  - A is (p*q) x (p*q), partitioned into p x p blocks of size q x q
#  - R(A) is p^2 x q^2 and must satisfy R(kron(B,C)) = vec(B) %*% t(vec(C))

kron <- function(A, B) kronecker(A, B)  # wrapper

# corrected rearrangement: row index r = (i2-1)*p + i1
rearrange_R <- function(A, p, q) {
  if (!is.matrix(A)) stop("A must be a matrix")
  if (nrow(A) != p*q || ncol(A) != p*q) {
    stop(sprintf("A must be (p*q) x (p*q). Got %d x %d", nrow(A), ncol(A)))
  }
  Rmat <- matrix(0, nrow = p * p, ncol = q * q)
  for (i1 in 1:p) {
    rows <- ((i1 - 1) * q + 1):(i1 * q)
    for (i2 in 1:p) {
      cols <- ((i2 - 1) * q + 1):(i2 * q)
      block_q_q <- A[rows, cols, drop = FALSE]   # q x q block
      # <-- fix: i2 contributes the slow index (column of B), i1 the fast index (row of B)
      Rrow <- (i2 - 1) * p + i1
      Rmat[Rrow, ] <- as.vector(block_q_q)       # column-major vec of the qxq block
    }
  }
  Rmat
}

# inverse vec (column-major)
unvec <- function(v, nr, nc) {
  if (length(v) != nr * nc) stop("length(v) != nr*nc")
  matrix(v, nrow = nr, ncol = nc)
}

# Main solver
best_kron_rank1 <- function(A, p, q, compute_full = FALSE) {
  if (!is.matrix(A)) stop("A must be a matrix.")
  if (nrow(A) != p*q || ncol(A) != p*q) stop("A must be (p*q)x(p*q).")
  
  Rm <- rearrange_R(A, p, q)
  s <- svd(Rm)  # need all singular values for separability
  sigma1 <- s$d[1]
  u1 <- s$u[, 1]
  v1 <- s$v[, 1]
  
  sqrt_sig <- sqrt(abs(sigma1))
  Bhat <- unvec(sqrt_sig * u1, p, p)
  Chat <- unvec(sqrt_sig * v1, q, q)
  Ahat <- kron(Bhat, Chat)
  residual_norm <- norm(A - Ahat, type = "F")
  
  # separability indices
  sep_index <- sigma1 / sum(s$d)
  sep_norm <- residual_norm / norm(A, "F")
  
  list(
    B = Bhat,
    C = Chat,
    sigma = sigma1,
    sep_index = sep_index,
    sep_norm = sep_norm,
    residual_norm = residual_norm,
    A_approx = if (compute_full) Ahat else NULL,
    Rmatrix = if (compute_full) Rm else NULL
  )
}

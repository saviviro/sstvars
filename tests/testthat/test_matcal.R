context("matcal")
library(sstvars)

a1 <- 1
a2 <- 1:4
a3 <- 1:7^2
A1 <- matrix(a1, nrow=1)
A2 <- matrix(a2, nrow=2, byrow=FALSE)
A3 <- matrix(a3, nrow=7, byrow=FALSE)


test_that("vec and unvec work correctly", {
  expect_equal(vec(A1), a1)
  expect_equal(vec(A2), a2)
  expect_equal(vec(A3), a3)
  expect_equal(unvec(d=1, a=a1), A1)
  expect_equal(unvec(d=2, a=a2), A2)
  expect_equal(unvec(d=7, a=a3), A3)
})

B1 <- matrix(1)
B2 <- matrix(c(1, 0.5, 0.5, 2), nrow=2, byrow=FALSE)
B3 <- matrix(c(1, 0.3, 0.2, 0.3, 2, 0.4, 0.2, 0.4, 3), nrow=3, byrow=FALSE)
b1 <- 1
b2 <- c(1, 0.5, 2)
b3 <- c(1, 0.3, 0.2, 2, 0.4, 3)

test_that("vech and unvech work correctly", {
  expect_equal(vech(B1), b1)
  expect_equal(vech(B2), b2)
  expect_equal(vech(B3), b3)
  expect_equal(unvech(d=1, a=b1), B1)
  expect_equal(unvech(d=2, a=b2), B2)
  expect_equal(unvech(d=3, a=b3), B3)
})

B1 <- matrix(0:3, nrow=2, byrow=FALSE)
B2 <- matrix(c(1, 2, 0, 3, 0, 4, 5, 6, 7), nrow=3, byrow=FALSE)
B3 <- matrix(c(1, 0, 2, 3, 4, 0, 0, 5, 0, 6, 7, 8, 0, 9, 0, 0), nrow=4, byrow=FALSE)
B4 <- matrix(1:25, nrow=5, byrow=FALSE)

test_that("Bvec and unBvec work correctly", {
  expect_equal(Bvec(B1), 1:3)
  expect_equal(Bvec(B2), 1:7)
  expect_equal(Bvec(B3), 1:9)
  expect_equal(Bvec(B4), 1:25)

  expect_equal(unBvec(Bvector=Bvec(B1), d=2, B_constraints=B1), B1)
  expect_equal(unBvec(Bvector=Bvec(B2), d=3, B_constraints=B2), B2)
  expect_equal(unBvec(Bvector=Bvec(B3), d=4, B_constraints=B3), B3)
  expect_equal(unBvec(Bvector=Bvec(B4), d=5, B_constraints=B4), B4)
})

Omega1_2 <- matrix(c(0.93, -0.15, -0.15, 5.20), nrow=2, byrow=FALSE) # d=2
Omega2_2 <- matrix(c(5.88, 3.56, 3.56, 9.80), nrow=2, byrow=FALSE)

Omega1_3 <- matrix(c(1, 0.22, 0.33, 0.22, 2, 0.44, 0.33, 0.44, 3), nrow=3, byrow=FALSE)
Omega2_3 <- matrix(c(1.1, 0.222, 0.333, 0.222, 2.2, 0.444, 0.333, 0.444, 3.3), nrow=3, byrow=FALSE)

make_W <- function(x, d) matrix(x[1:(d^2)], nrow=d, ncol=d, byrow=FALSE)
make_Lambda <- function(x, d) diag(x[(d^2 + 1):length(x)])

test_that("diag_Omegas works correctly", {
  x2 <- diag_Omegas(Omega1_2, Omega2_2)
  W2 <- make_W(x2, d=2)
  Lambda2 <- make_Lambda(x2, d=2)
  expect_equal(tcrossprod(W2), Omega1_2, tol=1e-6)
  expect_equal(W2%*%tcrossprod(Lambda2, W2), Omega2_2, tol=1e-6)

  x3 <- diag_Omegas(Omega1_3, Omega2_3)
  W3 <- make_W(x3, d=3)
  Lambda3 <- make_Lambda(x3, d=3)
  expect_equal(tcrossprod(W3), Omega1_3, tol=1e-6)
  expect_equal(W3%*%tcrossprod(Lambda3, W3), Omega2_3, tol=1e-6)
})

get_Omega <- function(M, d, W, lambdas, which, W_and_lambdas=NULL) {
  if(!is.null(W_and_lambdas)) {
    W <- W_and_lambdas[1:d^2]
    lambdas <- W_and_lambdas[(d^2 + 1):length(W_and_lambdas)]
  }
  W <- matrix(W, nrow=d, ncol=d, byrow=FALSE)
  if(which == 1) {
    return(tcrossprod(W))
  }
  lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE)
  W%*%tcrossprod(diag(lambdas[, which - 1]), W)
}

# Md
W22 <- 1:4
lambdas22 <- 1:2
Omega22_1 <- get_Omega(M=2, d=2, W=W22, lambdas=lambdas22, which=1)
Omega22_2 <- get_Omega(M=2, d=2, W=W22, lambdas=lambdas22, which=2)

W23 <- 1:9
lambdas23 <- 1:3
Omega23_1 <- get_Omega(M=2, d=3, W=W23, lambdas=lambdas23, which=1)
Omega23_2 <- get_Omega(M=2, d=3, W=W23, lambdas=lambdas23, which=2)

W32 <- (1:4)/2
lambdas32 <- (1:4)/5
Omega32_1 <- get_Omega(M=3, d=2, W=W32, lambdas=lambdas32, which=1)
Omega32_2 <- get_Omega(M=3, d=2, W=W32, lambdas=lambdas32, which=2)
Omega32_3 <- get_Omega(M=3, d=2, W=W32, lambdas=lambdas32, which=3)

W33 <- (1:9)/2
lambdas33 <- (1:6)/4
Omega33_1 <- get_Omega(M=3, d=3, W=W33, lambdas=lambdas33, which=1)
Omega33_2 <- get_Omega(M=3, d=3, W=W33, lambdas=lambdas33, which=2)
Omega33_3 <- get_Omega(M=3, d=3, W=W33, lambdas=lambdas33, which=3)

test_that("redecompose_Omegas works correctly", {
  # M=2, d=2
  decomp22_1 <- redecompose_Omegas(M=2, d=2, W=W22, lambdas=lambdas22, perm=1:2)
  decomp22_2 <- redecompose_Omegas(M=2, d=2, W=W22, lambdas=lambdas22, perm=2:1)

  expect_equal(get_Omega(M=2, d=2, W_and_lambdas=decomp22_1, which=1), Omega22_1, tolerance=1e-6)
  expect_equal(get_Omega(M=2, d=2, W_and_lambdas=decomp22_1, which=2), Omega22_2, tolerance=1e-6)
  expect_equal(get_Omega(M=2, d=2, W_and_lambdas=decomp22_2, which=1), Omega22_2, tolerance=1e-6)
  expect_equal(get_Omega(M=2, d=2, W_and_lambdas=decomp22_2, which=2), Omega22_1, tolerance=1e-6)

  # M=2, d=3
  decomp23 <- redecompose_Omegas(M=2, d=3, W=W23, lambdas=lambdas23, perm=2:1)

  expect_equal(get_Omega(M=2, d=3, W_and_lambdas=decomp23, which=1), Omega23_2, tolerance=1e-6)
  expect_equal(get_Omega(M=2, d=3, W_and_lambdas=decomp23, which=2), Omega23_1, tolerance=1e-6)

  # M=3, d=2
  decomp32_1 <- redecompose_Omegas(M=3, d=2, W=W32, lambdas=lambdas32, perm=c(1, 3, 2))
  decomp32_2 <- redecompose_Omegas(M=3, d=2, W=W32, lambdas=lambdas32, perm=c(2, 3, 1))
  decomp32_3 <- redecompose_Omegas(M=3, d=2, W=W32, lambdas=lambdas32, perm=c(3, 2, 1))

  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_1, which=1), Omega32_1, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_1, which=2), Omega32_3, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_1, which=3), Omega32_2, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_2, which=1), Omega32_2, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_2, which=2), Omega32_3, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_2, which=3), Omega32_1, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_3, which=1), Omega32_3, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_3, which=2), Omega32_2, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=2, W_and_lambdas=decomp32_3, which=3), Omega32_1, tolerance=1e-6)

  # M=3, d=3
  decomp33_1 <- redecompose_Omegas(M=3, d=3, W=W33, lambdas=lambdas33, perm=c(3, 1, 2))
  decomp33_2 <- redecompose_Omegas(M=3, d=3, W=W33, lambdas=lambdas33, perm=c(1, 3, 2))
  decomp33_3 <- redecompose_Omegas(M=3, d=3, W=W33, lambdas=lambdas33, perm=c(2, 1, 3))
  decomp33_4 <- redecompose_Omegas(M=3, d=3, W=W33, lambdas=lambdas33, perm=c(1, 2, 3))

  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_1, which=1), Omega33_3, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_1, which=2), Omega33_1, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_1, which=3), Omega33_2, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_2, which=1), Omega33_1, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_2, which=2), Omega33_3, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_2, which=3), Omega33_2, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_3, which=1), Omega33_2, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_3, which=2), Omega33_1, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_3, which=3), Omega33_3, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_4, which=1), Omega33_1, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_4, which=2), Omega33_2, tolerance=1e-6)
  expect_equal(get_Omega(M=3, d=3, W_and_lambdas=decomp33_4, which=3), Omega33_3, tolerance=1e-6)
})


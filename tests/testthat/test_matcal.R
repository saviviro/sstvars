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

test_that("Wvec and unWvec work correctly", {
  expect_equal(Wvec(B1), 1:3)
  expect_equal(Wvec(B2), 1:7)
  expect_equal(Wvec(B3), 1:9)
  expect_equal(Wvec(B4), 1:25)

  expect_equal(unWvec(Wvector=Wvec(B1), d=2, B_constraints=B1), B1)
  expect_equal(unWvec(Wvector=Wvec(B2), d=3, B_constraints=B2), B2)
  expect_equal(unWvec(Wvector=Wvec(B3), d=4, B_constraints=B3), B3)
  expect_equal(unWvec(Wvector=Wvec(B4), d=5, B_constraints=B4), B4)
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


test_that("get_symmetric_sqrt works correctly", {
  x2 <- get_symmetric_sqrt(Omega1_2)
  W2 <- make_W(x2, d=2)
  expect_equal(tcrossprod(W2), Omega1_2, tol=1e-6)

  x3 <- get_symmetric_sqrt(Omega1_3)
  W3 <- make_W(x3, d=3)
  expect_equal(tcrossprod(W3), Omega1_3, tol=1e-6)
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



mat_power_test <- function(A, j) {
  if(j == 0) return(diag(nrow(A)))
  Reduce('%*%', replicate(j, A, simplify=FALSE))
}

B5 <- cbind(rbind(B2, Omega1_3), rbind(Omega23_1, Omega1_3))
W3 <- matrix(W32, nrow=2, ncol=2)
W4 <- matrix(W33, nrow=3, ncol=3)


test_that("mat_power works correctly", {
  expect_equal(mat_power(A1, 0), A1, tol=1e-10)
  expect_equal(mat_power(A1, 1), A1, tol=1e-10)
  expect_equal(mat_power(A1, 2), A1%*%A1, tol=1e-10)

  expect_equal(mat_power(A2, 0), diag(nrow(A2)), tol=1e-10)
  expect_equal(mat_power(A2, 1), A2, tol=1e-10)
  expect_equal(mat_power(A2, 2), A2%*%A2, tol=1e-10)
  expect_equal(mat_power(A2, 3), A2%*%A2%*%A2, tol=1e-10)
  expect_equal(mat_power(A2, 13), mat_power_test(A2, 13), tol=1e-10)
  expect_equal(mat_power(A2, 51), mat_power_test(A2, 51), tol=1e-10)

  expect_equal(mat_power(A3, 0), diag(nrow(A3)), tol=1e-10)
  expect_equal(mat_power(A3, 1), A3, tol=1e-10)
  expect_equal(mat_power(A3, 2), A3%*%A3, tol=1e-10)
  expect_equal(mat_power(A3, 3), A3%*%A3%*%A3, tol=1e-10)
  expect_equal(mat_power(A3, 8), mat_power_test(A3, 8), tol=1e-10)
  expect_equal(mat_power(A3, 27), mat_power_test(A3, 27), tol=1e-10)

  expect_equal(mat_power(B3, 3), B3%*%B3%*%B3, tol=1e-10)
  expect_equal(mat_power(B3, 4), mat_power_test(B3, 4), tol=1e-10)
  expect_equal(mat_power(B3, 43), mat_power_test(B3, 43), tol=1e-10)

  expect_equal(mat_power(B4, 0), diag(nrow(B4)), tol=1e-10)
  expect_equal(mat_power(B4, 1), B4, tol=1e-10)
  expect_equal(mat_power(B4, 2), B4%*%B4, tol=1e-10)
  expect_equal(mat_power(B4, 3), B4%*%B4%*%B4, tol=1e-10)
  expect_equal(mat_power(B4, 4), B4%*%B4%*%B4%*%B4, tol=1e-10)
  expect_equal(mat_power(B4, 10), mat_power_test(B4, 10), tol=1e-10)
  expect_equal(mat_power(B4, 21), mat_power_test(B4, 21), tol=1e-10)

  expect_equal(mat_power(B5, 0), diag(nrow(B5)), tol=1e-10)
  expect_equal(mat_power(B5, 1), B5, tol=1e-10)
  expect_equal(mat_power(B5, 2), B5%*%B5, tol=1e-10)
  expect_equal(mat_power(B5, 3), B5%*%B5%*%B5, tol=1e-10)
  expect_equal(mat_power(B5, 4), B5%*%B5%*%B5%*%B5, tol=1e-10)
  expect_equal(mat_power(B5, 10), mat_power_test(B5, 10), tol=1e-10)
  expect_equal(mat_power(B5, 21), mat_power_test(B5, 21), tol=1e-10)

  expect_equal(mat_power(W3, 0), diag(nrow(W3)), tol=1e-10)
  expect_equal(mat_power(W3, 1), W3, tol=1e-10)
  expect_equal(mat_power(W3, 2), W3%*%W3, tol=1e-10)
  expect_equal(mat_power(W3, 3), W3%*%W3%*%W3, tol=1e-10)
  expect_equal(mat_power(W3, 5), mat_power_test(W3, 5), tol=1e-10)
  expect_equal(mat_power(W3, 24), mat_power_test(W3, 24), tol=1e-10)

  expect_equal(mat_power(W4, 0), diag(nrow(W4)), tol=1e-10)
  expect_equal(mat_power(W4, 1), W4, tol=1e-10)
  expect_equal(mat_power(W4, 2), W4%*%W4, tol=1e-10)
  expect_equal(mat_power(W4, 3), W4%*%W4%*%W4, tol=1e-10)
  expect_equal(mat_power(W4, 4), W4%*%W4%*%W4%*%W4, tol=1e-10)
  expect_equal(mat_power(W4, 11), mat_power_test(W4, 11), tol=1e-10)
  expect_equal(mat_power(W4, 20), mat_power_test(W4, 20), tol=1e-10)
})

test_that("create_J_matrix works", {

  # Test with d = 2, p = 3
  J <- create_J_matrix(2, 3)

  # Check dimensions
  expect_equal(dim(J), c(2, 6))

  # Check the first d x d block is identity matrix
  expect_equal(J[ , 1:2], diag(2))

  # Check remaining blocks are zeros
  expect_equal(J[ , 3:6], matrix(0, nrow = 2, ncol = 4))

  # Test with d = 3, p = 1 (edge case where p = 1)
  J <- create_J_matrix(3, 1)

  # Check dimensions
  expect_equal(dim(J), c(3, 3))

  # Check it's an identity matrix
  expect_equal(J, diag(3))

  # Test with d = 100, p = 50 (large case)
  J <- create_J_matrix(100, 50)

  # Check dimensions
  expect_equal(dim(J), c(100, 5000))

  # Check the first d x d block is identity matrix
  expect_equal(J[ , 1:100], diag(100))

  # Check remaining blocks are zeros
  expect_equal(J[ , 101:5000], matrix(0, nrow = 100, ncol = 4900))

  # Test with d = 5, p = 9 (moderate size case)
  J <- create_J_matrix(5, 9)

  # Check dimensions
  expect_equal(dim(J), c(5, 45))

  # Check the first K x K block is identity matrix
  expect_equal(J[ , 1:5], diag(5))

  # Check remaining blocks are zeros
  expect_equal(J[ , 6:45], matrix(0, nrow = 5, ncol = 40))
})


# Test that the function creates a matrix with the correct dimensions
test_that("create_Fi_matrix returns matrix with correct dimensions", {
  for(T_obs in c(5, 10, 20, 213)) {
    for(i in 0:(T_obs - 1)) {
      expect_equal(dim(create_Fi_matrix(i=i, T_obs=T_obs)), c(T_obs, T_obs))
    }
  }
})

# Test that the sub-diagonal is correctly filled with 1's
test_that("Sub-diagonal of Fi is correctly filled with 1's", {
  for(T_obs in c(5, 10, 16, 19)) {
    for(i in 1:(T_obs - 1)) {
      Fi <- create_Fi_matrix(i=i, T_obs=T_obs)
      expected_Fi <- matrix(data=0, nrow=T_obs, ncol=T_obs)
      expected_Fi[cbind((i+1):T_obs, 1:(T_obs-i))] <- 1
      expect_equal(Fi, expected_Fi)
    }
  }
})

# Test that the function returns an identity matrix for i=0
test_that("create_Fi_matrix returns identity matrix for i=0", {
  for(T_obs in c(1, 5, 10, 13)) {
    Fi <- create_Fi_matrix(i=0, T_obs=T_obs)
    expected_Fi <- diag(x=1, nrow=T_obs)
    expect_equal(Fi, expected_Fi)
  }
})

# Test that the function returns a zero matrix for T_obs equal to i
test_that("create_Fi_matrix returns zero matrix for T_obs equal to i", {
  for(T_obs in c(1, 5, 10, 13, 20)) {
    Fi <- create_Fi_matrix(i=T_obs, T_obs=T_obs)
    expected_Fi <- matrix(data=0, nrow=T_obs, ncol=T_obs)
    expect_equal(Fi, expected_Fi)
  }
})


test_that("order_B works correctly", {
  B <- matrix(c(-1, 2, 3, -4), nrow=2)
  expected <- matrix(c(3, -4, 1, -2), nrow=2)
  expect_equal(order_B(B), expected)

  B <- matrix(c(-2, -3, 4, 5), nrow=2)
  expected <- matrix(c(4, 5, 2, 3), nrow=2)
  expect_equal(order_B(B), expected)
  expect_equal(order_B(expected), expected)

  B <- matrix(c(-10, 20, -30, 40, -50, 60, -70, 80, 90), nrow=3)
  expected <- matrix(c(70, -80, -90, 40, -50, 60, 10, -20, 30), nrow=3)
  expect_equal(order_B(B), expected)

  B <- matrix(c(-50, 20, -30, 40, -50, 60, -70, 80, 90), nrow=3)
  expected <- matrix(c(70, -80, -90, 50, -20, 30, 40, -50, 60), nrow=3)
  expect_equal(order_B(B), expected)

  B <- matrix(c(-10, 20, -30, 40, -50, 60, 20, 80, 90), nrow=3)
  expected <- matrix(c(40, -50, 60, 20, 80, 90, 10, -20, 30), nrow=3)
  expect_equal(order_B(B), expected)
})




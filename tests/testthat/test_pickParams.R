context("pickParams")
library(sstvars)

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(0.65, 0.7)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)

theta_112relg <- c(phi10_112, vec(A11_112), vech(Omega1_112))

# p=2, M=1, d=2
phi10_212 <- c(0.53, 0.03)
A11_212 <- matrix(c(0.23, 0.02, -0.17, 0.66), nrow=2, byrow=FALSE)
A12_212 <- matrix(c(0.18, 0.02, 0.04, 0.26), nrow=2, byrow=FALSE)
Omega1_212 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

theta_212relg <- c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212))

# p=1, M=2, d=2
phi10_122 <- c(0.55, 0.11)
A11_122 <- matrix(c(0.34, 0.05, -0.01, 0.72), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

phi20_122 <- c(0.17, 0.25)
A21_122 <- A11_122
Omega2_122 <- matrix(c(0.50, -0.01, -0.01, 0.20), nrow=2, byrow=FALSE)

alpha1_122 <- 0.60
theta_122relg <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)

# p=2, M=2, d=2
phi10_222 <- c(0.36, 0.12)
A11_222 <- matrix(c(0.22, 0.06, -0.15, 0.39), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(0.41, -0.01, 0.08, 0.3), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.21, 0.01, 0.01, 0.03), nrow=2, byrow=FALSE)

phi20_222 <- c(0.48, 0.07)
A21_222 <- matrix(c(0.22, 0.02, -0.12, 0.72), nrow=2, byrow=FALSE)
A22_222 <- matrix(c(0.09, 0.03, 0.04, 0.19), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(1.10, 0.01, 0.01, 0.11), nrow=2, byrow=FALSE)

alpha1_222 <- 0.37
theta_222relg <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                   vech(Omega1_222), vech(Omega2_222), alpha1_222)

# p=1, M=3, d=2
phi10_132 <- phi10_122
phi20_132 <- phi20_122
phi30_132 <- c(12, 13)

A11_132 <- A11_122
A21_132 <- A21_122
A32_132 <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow=2)
Omega1_132 <- Omega1_122
Omega2_132 <- Omega2_122
Omega3_132 <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
alpha1_132 <- 0.5
alpha2_132 <- 0.3

theta_132relg <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A21_132),
                   vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

# p=1, M=1, d=3
phi10_113 <- c(1, 2, 3)
A11_113 <- matrix(c(0.1, 0.02, 0.12, 0.3, 0.21, 0.11, 0.05, 0.03, 0.09), nrow=3)
Omega1_113 <- matrix(c(c(1, 0.2, 0.3, 0.2, 2, 0.4, 0.3, 0.4, 3)), nrow=3)

theta_113relg <- c(phi10_113, vec(A11_113), vech(Omega1_113))

test_that("pick_phi0 work correctly", {
  expect_equal(pick_phi0(M=1, d=2, params=theta_112relg), as.matrix(phi10_112))
  expect_equal(pick_phi0(M=1, d=2, params=theta_212relg), as.matrix(phi10_212))
  expect_equal(pick_phi0(M=2, d=2, params=theta_122relg)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122relg)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222relg)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222relg)[,2], phi20_222)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relg)[,1], phi10_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relg)[,2], phi20_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relg)[,3], phi30_132)
  expect_equal(pick_phi0(M=1, d=3, params=theta_113relg), as.matrix(phi10_113))
})


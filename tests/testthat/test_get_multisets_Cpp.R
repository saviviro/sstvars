context("get_permanent_Cpp")
library(sstvars)

mat1 <- matrix(c(1, 1,
                 1, 2,
                 2, 2), nrow=3, byrow=TRUE)
mat2 <- matrix(c(1, 1, 1,
                 1, 1, 2,
                 1, 2, 2,
                 2, 2, 2), nrow=4, byrow=TRUE)
mat3 <- matrix(c(1, 1,
                 1, 2,
                 1, 3,
                 2, 2,
                 2, 3,
                 3, 3), nrow=6, byrow=TRUE)
mat4 <- matrix(c(1, 1, 1,
                 1, 1, 2,
                 1, 1, 3,
                 1, 2, 2,
                 1, 2, 3,
                 1, 3, 3,
                 2, 2, 2,
                 2, 2, 3,
                 2, 3, 3,
                 3, 3, 3), nrow=10, byrow=TRUE)
mat5 <- matrix(c(1, 1, 1, 1,
                 1, 1, 1, 2,
                 1, 1, 1, 3,
                 1, 1, 2, 2,
                 1, 1, 2, 3,
                 1, 1, 3, 3,
                 1, 2, 2, 2,
                 1, 2, 2, 3,
                 1, 2, 3, 3,
                 1, 3, 3, 3,
                 2, 2, 2, 2,
                 2, 2, 2, 3,
                 2, 2, 3, 3,
                 2, 3, 3, 3,
                 3, 3, 3, 3), nrow=15, byrow=TRUE)

test_that("get_multisets_Cpp work correctly", {
  expect_equal(get_multisets_Cpp(n=2, d=1, N=choose(2 + 1 - 1, 1)), matrix(1:2, nrow=2), tolerance=1e-3)
  expect_equal(get_multisets_Cpp(n=2, d=2, N=choose(2 + 2 - 1, 2)), mat1, tolerance=1e-3)
  expect_equal(get_multisets_Cpp(n=2, d=3, N=choose(2 + 3 - 1, 3)), mat2, tolerance=1e-3)
  expect_equal(get_multisets_Cpp(n=3, d=2, N=choose(3 + 2 - 1, 2)), mat3, tolerance=1e-3)
  expect_equal(get_multisets_Cpp(n=3, d=3, N=choose(3 + 3 - 1, 3)), mat4, tolerance=1e-3)
  expect_equal(get_multisets_Cpp(n=3, d=4, N=choose(3 + 4 - 1, 4)), mat5, tolerance=1e-3)
})



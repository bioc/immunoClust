
test_that("bhattacharyya", {
  
  S <- diag(0.5,2)
  M1 <- c(0,0)
  M2 <- c(1,1)
  #ld <- det(S) = 0.25
  d <- bhattacharyya.dist(M1,S,M2,S)
  c <- bhattacharyya.coeff(M1,S,M2,S)
  p <- bhattacharyya.prob(M1,S,M2,S)
  
  expect_equal(d, 0.5)
  expect_equal(-log(c), d)
  expect_equal(log(p)+0.125*log(0.25), -d)
  
})




test_that("transformation", {
  
  data(dat.fcs)
  data(dat.exp)
  dat.trans <- trans.ApplyToData(dat.exp[[1]], dat.fcs)
  expect_is(dat.trans, "flowFrame")
  r <- cell.Classify(dat.exp[[1]], dat.trans)
  expect_is(r, "immunoClust")
  t <- trans.FitToData(r, dat.fcs)
  expect_length(t,14)
})



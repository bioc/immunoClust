
test_that("transformation", {
  
 data(dat.fcs)
 data(dat.exp)
 dat.trans <- trans.ApplyToData(dat.exp[[1]], dat.fcs)
 
 r <- cell.Classify(dat.exp[[1]], dat.trans)
 t <- trans.FitToData(r, dat.fcs)
  
})



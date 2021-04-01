
test_that("cell.EM", {

    data(dat.fcs)
    data(dat.exp)
    ## cell.clustering result for dat.fcs
    r <- dat.exp[[1]]
    ## apply model parameter to all (unfiltered) events
    dat.trans <- trans.ApplyToData(r, dat.fcs)
    r2 <- cell.EM(dat.trans, parameters(r), K=ncls(r),
                w=weights(r),m=mu(r),s=sigma(r))

    expect_is(r2, "immunoClust")

})

test_that("cell.FitModel", {

    data(dat.fcs)
    data(dat.exp)
    r1 <- dat.exp[[1]]
    dat.trans <- trans.ApplyToData(r1, dat.fcs)
    r2 <- cell.FitModel(r1, dat.trans)

    expect_is(r2, "immunoClust")
})

test_that("cell.ClustData", {

    data(dat.fcs)
    res <- cell.ClustData(dat.fcs, 5, parameters=c("FSC-A", "SSC-A"))

    expect_is(res, "immunoClust")
})

test_that("cell.hclust", {

    data(dat.fcs)
    inc <- sample(1:nrow(dat.fcs), 50)
    result <- cell.hclust(exprs(dat.fcs)[inc,])
    expect_length(result, 98)
})


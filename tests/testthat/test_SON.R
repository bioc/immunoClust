test_that("meta.SON.clustering", {
    data(dat.meta)
    res <- meta.SON.clustering(dat.meta, cycles=2)
    expect_is(res, "immunoMeta")
})

test_that("meta.SON.combineClustering", {
    data(dat.exp)
    data(dat.meta)
    res <- meta.SON.combineClustering(dat.meta, dat.exp[[1]], SON.cycles=2)
    expect_is(res, "immunoClust")
})

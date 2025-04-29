
progs = c("pedprobr", "merlin", "Familias", "ElstonStewart")
p = 0.9; q = 1-p

test_that("all programs work", {
  x = nuclearPed(1) |> addMarker()
  res = compare(x, verbose=F)
  expect_true(all_agree(res, answer=1))
})

test_that("all programs agree in inbred example", {
  x = halfCousinPed(0, child = T) |>
    addMarker('4' = "2/2", '6' = "2/2", alleles=1:2, afreq=c(p,q))
  res = compare(x, verbose=F)
  expect_true(all_agree(res))
})

test_that("all programs agree in simple X-linked example", {
  x = nuclearPed(1) |>
    addMarker('3' = 2, alleles=1:2, afreq=c(p,q), chrom=23)
  res = compare(x, verbose=F)
  expect_true(all_agree(res, answer=q))
})

test_that("all programs agree in inbred X-linked", {
  x = halfCousinPed(0, child=T) |>
    addMarker('6' = 2, alleles=1:2, afreq=c(p,q), chrom=23)
  res = compare(x, verbose=F)
  expect_true(all_agree(res, q))
})

test_that("errors are caught", {
  skip("Weird")
  x = quadHalfFirstCousins() |>
    addMarker('10' = "1/2")
  expect_error(likelihood_pedprobr(x, verbose=F),
               "This pedigree requires founders as loop breakers")
  res = compare(x, verbose=F)
  expect_equal(res$program, c("Familias","ElstonStewart","merlin"))
  expect_true(all_agree(res))
})

test_that("mutation modelling works", {
  x = nuclearPed() |>
    addMarker("1" = "1/1", "3" = "2/2", alleles = 1:5, mutmod = "eq", rate = 0.1)
  res = compare(x, verbose=F)
  expect_true(all_agree(res))
})


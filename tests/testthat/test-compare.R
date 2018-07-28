context("compare")

progs = c("pedprobr", "paramlink", "merlin", "Familias", "ElstonStewart")

test_that("all programs work", {
  x = nuclearPed(1)
  m = marker(x, alleles=1:2)
  res = compare(x, m, verbose=F)
  expect_true(all_agree(res, answer=1))
})

test_that("all programs agree in inbred example", {
  x = halfCousinsPed(0, child=T)
  p = 0.9; q = 1-p
  m = marker(x, '6' = 2, alleles=1:2, afreq=c(p,q))
  res = compare(x, m, verbose=F)
  expect_true(all_agree(res, answer=1/8*q + 7/8*q^2))
})

test_that("all programs agree in simple X-linked example", {
  x = nuclearPed(1)
  p = 0.9; q = 1-p
  m = marker(x, '3' = 2, alleles=1:2, afreq=c(p,q), chrom=23)
  res = compare(x, m, verbose=F)
  expect_true(all_agree(res, answer=q))
})

test_that("all programs agree in inbred X-linked", {
  x = halfCousinsPed(0, child=T)
  p = 0.9; q = 1-p
  m = marker(x, '6' = 2, alleles=1:2, afreq=c(p,q), chrom=23)
  res = compare(x, m, verbose=F)
  expect_true(all_agree(res, q))
})

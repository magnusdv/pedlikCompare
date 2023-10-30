#' Check if all programs agree
#'
#' Analyses the output table of [compare()] and checks if all programs agree.
#' MERLIN is treated separately, since it only gives an approximate answer. (It
#' reports the log-likelihood rounded to 3 decimals.)
#'
#' @param df A `data.frame` (or `tibble`) produced by `compare()`
#' @param answer The correct likelihood. If not provided, the first value in the
#'   `likelihood` column is used.
#'
#' @return TRUE or FALSE
#'
#' @examples
#'
#' x = nuclearPed(1) |> addMarker(`1` = "1/2")
#' res = compare(x)
#' stopifnot(all_agree(res, answer = 0.5))
#'
#' @export
all_agree = function(df, answer=NULL) {
  if(is.null(answer))
    answer = df$lnlik[1]
  else
    answer = log(answer)

  df_without_merlin = df[df$program != "merlin", ]
  test1 = isTRUE(all.equal(df_without_merlin$lnlik,
                           rep(answer, nrow(df_without_merlin))))

  test2 = TRUE

  # Test MERLIN
  if("merlin" %in% df$program) {
    merlin_lnlik = df$lnlik[match("merlin", df$program)]
    deci = decims(merlin_lnlik)
    test2 = isTRUE(all.equal(merlin_lnlik, round(answer, deci)))
  }

  test1 && test2
}

# Hack to find the number of decimals in a number
decims = function(x) {
  tol = .Machine$double.eps^0.5

  shift = x*10^(0:20)
  min(which(abs(shift - round(shift)) < tol)) - 1
}

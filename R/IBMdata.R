#' Daily IBM log-returns.
#'
#' A dataset containing the daily IBM log-returns \eqn{r_t=log(s_t/s_{t-1})},
#' where \eqn{s_t} denotes the stock price (closing to closing) at day \eqn{t}
#' from Jan. 3rd, 2006 - Dec. 31, 2018 with overall series length \eqn{T=3271}.
#'
#' @format A data frame with 3271 rows and 1 variable:
#' \describe{
#'   \item{r_t}{daily log-return}
#' }
#' @usage data(IBMdata)
#' @source \url{https://finance.yahoo.com/}
"IBMdata"

#' Enumerates potential matches
#' 
#' Given a variable \code{Treat}, \code{matching_group} enumerates the individuals
#' towards each observation can be matched against. If \code{exact} is specified,
#' then it performes exact match (see examples).
#' 
#' @template mat
#' @templateVar exact 1
#' @templateVar Treat 1
#' @return A list of length \eqn{n} specifying for each i to which individuals
#' each j != i will be matched.
#' 
#' @details In the most common case \code{Treat} will be either of \code{0} or
#' \code{1}, in which case individuals will be matched against their opposite
#' group. If \code{Treat == -1}, then the group constraint will be void.
#' 
#' @name matching_group
#' @examples
#' # Asigning matching groups for lalonde --------------------------------------
#' 
#' data(lalonde, package = "MatchIt")
#' dat <- lalonde
#' 
#' # Match (no exact)
#' m_noexact <- matching_group(dat$treat)
#' 
#' # How many matches?
#' table(sapply(m_noexact, length))
#' table(dat$treat)
#' 
#' # What if we ask exact match on black ---------------------------------------
#' 
#' m_exact <- matching_group(dat$treat, dat$black)
#' 
#' table(sapply(m_exact, length))
#' with(dat, table(treat, black))
#' @export
matching_group <- function(Treat, exact = NULL) {
  
  # Checking exact
  if (!length(exact))
    exact <- cbind(rep(1, length(Treat)))
  else 
    exact <- as.matrix(exact)
  
  # Checking treat
  Treat <- cbind(Treat)
  
  matching_group_cpp(Treat, exact, 1)
}
#' Identify duplicated rows and group them keeping track of memberships
#' 
#' A lot of times observations may have duplicates within a same group.
#' For example, binary covariates have a high probability to be repeated.
#' Internally, [blop] controls for this and when performing the matching
#' takes into account duplicates so that the weights are equally distributed
#' accross duplicates.
#' 
#' @param x Either a Numeric matrix or a data-frame.
#' @inheritParams blop
#' @param ids Integer vector with ids to be mapped to the set (internal use).
#' 
#' @return A nested list with the following elements:
#' - `groups` Lists of integers indicating the ids (row position) of
#' the group members.
#' - `selected` Integer vector of length `length(groups)` with a proposed
#' set of ids to be used for the matching.
#' 
#' If `Treat` is passed, then it will return a list of length 3 with nested
#' lists as described above. 
#' 
#' @examples 
#' data(nsw)
#' vars <- c("black", "hispanic", "nodegree", "re75")
#' str(group_duplicates(nsw[1:10, vars]))
#' str(group_duplicates(nsw[290:310, vars], Treat = nsw$treat[290:310]))
#' @export
group_duplicates <- function(x, Treat = NULL, ids = NULL) UseMethod("group_duplicates")

#' @export
#' @rdname group_duplicates
group_duplicates.data.frame <- function(x, Treat = NULL, ids = NULL) {
  
  # Identifying string variables
  for (v in 1L:ncol(x)) {
    
    # Case of character variables
    if (inherits(x[[v]], "character"))
      x[[v]] <- as.factor(x[[v]])
    
    # Coercing into numeric
    x[[v]] <- as.numeric(x[[v]])
  }
  
  group_duplicates(as.matrix(x), Treat, ids)
  
  
}

#' @rdname group_duplicates
group_duplicates.matrix <- function(x, Treat = NULL, ids = NULL) {
  
  # Checking types
  if (typeof(x) == "logical")
    x[] <- as.integer(x)
  else if (!any(c("integer", "double") %in% typeof(x)))
    stop("-x- cannot be coerced into a numeric matrix.", call. = FALSE)
  
  # Checking NAs
  if (any(is.na(x)))
    stop("-x- has one or more NA values.", call. = FALSE)
  
  if (is.null(ids))
    ids <- 1L:nrow(x)
  
  # Checking Treat
  if (!is.null(Treat)) {
    
    # Same length?
    if (length(Treat) != nrow(x))
      stop("-Treat- must be of length nrow(x).", call. = FALSE)
    
    ttreat <- table(Treat)
    values <- -1:1
    if (!all(names(ttreat) %in% as.integer(values)))
      stop("-Treat- must have values in -1, 0, or 1.", call. = FALSE)
    
    ans <- structure(vector("list", 3), names = values)
    for (v in values) {
      
      # Checking the subset
      idx <- which(Treat == v)
      if (length(idx) == 0)
        next
      
      # Doing the matching within the group
      ans[[as.character(v)]] <- group_duplicates(
        x     = x[idx, , drop = FALSE],
        Treat = NULL,
        ids   = ids[idx]
        )
      
    }
    
    return(ans)
    
  }
  
  groups   <- vector("list", nrow(x))
  selected <- integer(0)
  xt       <- t(x)
  toskip   <- NULL
  count    <- 1L
  
  for (i in 1:nrow(x)) {
    
    if (ids[i] %in% toskip)
      next
    
    selected <- c(selected, ids[i])
    
    # Comparing
    m <- ids[which(colSums((xt - x[i,])^2) == 0)]
    m <- setdiff(m, c(ids[i], toskip))
    
    if (length(m) > 0) {
      
      groups[[count]] <- c(groups[[count]], m, ids[i])
      toskip <- c(toskip, m, ids[i])
      
    } else {
      
      groups[[count]] <- ids[i]
      toskip <- c(toskip, ids[i])
      
    }
    
    # Incrementing the count
    count <- count + 1L
  }
  
  list(
    groups   = groups[1:length(selected)],
    selected = selected
  )
}
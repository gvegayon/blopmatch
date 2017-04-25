#' \code{}
#' <%= ifelse(exists("X"), "@param X Numeric matrix of size \\eqn{n\\times K}{n*K}.", "") %>
#' <%= ifelse(exists("W"), "@param W Numeric matrix of size \\eqn{K\\times K}{K*K}.", "") %> 
#' <%= ifelse(exists("Treat"), "@param Treat Integer vector of size \\eqn{n}. Group indicator. Can be either -1 (no group), 0 (control) or 1 (treatment).", "") %> 
#' <%= ifelse(exists("exact"), "@param exact Numeric matrix with \\eqn{n} rows. When specified, exact matching is performed on those variables (see \\code{\\link{matching_group}}).", "")%>


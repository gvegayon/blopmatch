#' \code{}
#' <%= ifelse(exists("X"), "@param X Numeric matrix of size \\eqn{n\\times K}{n*K}.", "") %>
#' <%= ifelse(exists("y"), "@param y Numeric vector of length \\eqn{n}. Dependent variable.", "") %>
#' <%= ifelse(exists("W"), "@param W Numeric matrix of size \\eqn{K\\times K}{K*K}. When \\code{W} is equal to \\code{solve(var(X))}, it is the Mahalanobis norm.", "") %> 
#' <%= ifelse(exists("p"), "@param p Numeric scalar. See \\code{\\link{weighted_norm}}.", "") %>
#' <%= ifelse(exists("Treat"), "@param Treat Integer vector of size \\eqn{n}. Group indicator. Can be either -1 (no group), 0 (control) or 1 (treatment).", "") %> 
#' <%= ifelse(exists("exact"), "@param exact Numeric matrix with \\eqn{n} rows. When specified, exact matching is performed on those variables (see \\code{\\link{matching_group}}).", "")%>


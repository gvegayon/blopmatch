

# 1.  Clean the data: Leave 1 replicate per set.
#     a.  vector of classes
#     b.  

clean <- function(Treat, x) {
  
  # Pulling indices
  is1 <- Treat == 1
  is0 <- which(!is1)
  is1 <- which(is1)
  
  # Finding matches
  g0 <- matching_group(rep(-1, length(is0)), x[is0, , drop=FALSE])
  g1 <- matching_group(rep(-1, length(is1)), x[is1, , drop=FALSE])
  
  
  
}
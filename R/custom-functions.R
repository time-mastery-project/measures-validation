normalize = function(x, naIgnore=T){
  q = rep(NA,length(x))
  if(naIgnore==F & sum(is.na(x))>0) return(q)
  p = rank(x[!is.na(x)])/(length(x[!is.na(x)])+1)
  q[!is.na(x)] = qnorm(p)
  return(q)
}

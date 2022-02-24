pairwise_ks_test <- function(value, group, n_min = 50, warning = 0, alternative = "two.sided" ){
  
  lev <- unique(group)
  
  lst <- lapply( seq_along(lev), function(i) value[group == lev[i]] )
  names(lst)<-lev
  
  if (sum(lengths(lst)< n_min)) {
    lst <- lst [-which(lengths(lst)< n_min)]}
  
  f <- function(x, y){ 
    w <- getOption("warn") 
    options(warn = warning)  # ignore warnings 
    p <- ks.test(x, y, alternative = alternative, exact = 
                   F)$p.value 
    options(warn = w) 
    return(p) 
  } 
  
  res <- lapply(lst, function(x) lapply(lst, function(y) f(x, y))) 
  
  res<-unlist(res)
  res <- matrix(res, nrow = length(lst), ncol = length(lst), byrow = T)
  row.names(res) <- colnames(res) <- names(lst)
  cat("Pairwise Kolmogorov-Smirnov Test p-value Matrix","\n","\n")
  return(res)
}

prune <- function(pval,method, ...){
# Each method must output the number removed!
if(method=="partha"){
  percent_prune <- min(.2, ((pval - constant*control_limit)/(1-constant*control_limit))^2)
  percent_keep <- 1 - percent_prune
  num_keep <- ceiling(dim(data)[2] * percent_keep)
  num_remove <- dim(data)[2]-num_keep
}else if(method=="constant"){
  percent_prune <- constant # Note that constant in this context is quite different
  percent_keep <- 1 - percent_prune
  num_keep <- ceiling(dim(data)[2] * percent_keep)
  num_remove <- dim(data)[2]-num_keep
}else if(method=="ewma"){
  percent_prune <- .031  # hard coded to avoid doubling parameters right now
  percent_keep <- 1 - percent_prune
  num_keep <- ceiling(dim(data)[2] * percent_keep)
  num_remove <- dim(data)[2]-num_keep
}else if(method=="lookback"){
  # constant in this context is the max window size
  W = dim(data)[2] # note: can replace these calls to dim(data)[2] w/ num_active_batches
  if(pval < pval.small){
    num_remove = 2
  }else{
    num_remove = 0
  }
  if(W - num_remove > constant) num_remove = 1
  return(num_remove)
}else{
  stop("invalid method invoked")
}
}
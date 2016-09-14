sp.test.unif <- function(data,
                       levels,
                       usecrit = 0,
                       critvalL = 0,
                       critvalR = 0)
{
  res <- SP_Test_Uniformity(data)
  decisions <- rep(0, length(levels))
  for (i in 1:length(levels)) {
    if (usecrit == 0) {
      decisions[i] <- if (res$p.value < levels[i])
        1
      else
        0
    } else {
      decisions[i] <- if (res$statistic > critvalR[i])
        1
      else
        0
    }
  }
  return(
    list(
      statistic = res$statistic,
      pvalue = res$p.value,
      decision = decisions,
      alter = 3,
      stat.pars = NULL,
      pvalcomp = 0L,
      nbparstat = 0
    )
  )
}
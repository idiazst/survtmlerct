print.adjusted_rmst <- function(x){
  cat('Estimate for control arm (A == 0): ', x$rmst['rmst0'], '\n')
  cat('Estimate for treatment arm (A == 1): ', x$rmst['rmst1'], '\n') 
  cat('Standard error for the difference: ', x$std.error.diff, '\n') 
}

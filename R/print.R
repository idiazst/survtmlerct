print.adjusted_rmst <- function(x){
  cat('RMST for control arm (A == 0): ', x$rmst['rmst0'], '\n')
  cat('RMST for treatment arm (A == 1): ', x$rmst['rmst1'], '\n') 
  cat('Standard error for the difference: ', x$std.error.diff, '\n') 
}
print.ipw_rmst <- function(x){
  cat('RMST for control arm (A == 0): ', x['rmst0'], '\n')
  cat('RMST for treatment arm (A == 1): ', x'rmst1'], '\n') 
  cat('Standard errors not available for IPW\n') 
}
print.adjusted_prob <- function(x){
  cat('Survival probability for control arm (A == 0): ', x$prob['prob0'], '\n')
  cat('Survival probability for treatment arm (A == 1): ', x$prob['prob1'], '\n') 
  cat('Standard error for the difference: ', x$std.error.diff, '\n') 
}


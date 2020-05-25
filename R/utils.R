#' Bound a vector of probabilities in [r, 1 - r] for r > 0
#'
#' @param x Vector of probabilities
#' @param r Bound
#'
#' @return Bounded vector of probabilities
bound01 <- function(x, r = 1e-10){
    xx <- x
    xx[x < r] <- r
    xx[x > 1-r] <- 1-r
    return(as.numeric(xx))
}

#' Bound a vector of probabilities in [r, 1] for r > 0
#'
#' @param x Vector of probabilities
#' @param r Bound
#'
#' @return Bounded vector of probabilities
bound <- function(x, r = 0.001){
    xx <- x
    xx[x < r] <- r
    return(as.numeric(xx))
}

#' Transform a survival dataset from short to long form
#'
#' @param df A data frame containing columns T (time to event) and C (censoring indicator)
#' @param freq.time If time T must be ccategorized into coarser intervals (e.g., use freq.time = 7 if T is in days and you want the output in weeks)
#'
#' @return A long form data set
transformData <- function(df, freq.time){

    ## Transform a dataset from the short to the long form.
    if(freq.time > 1) df$T <- df$T %/% freq.time + 1

    n <- dim(df)[1]
    K <- max(df$T)

    m <- rep(1:K, n)
    Lm <- Rm <- rep(NA, n*K)
    Im <- Jm <- 1*(m == 1)

    ## Note: R and J are lagged with respect to definitions in the paper.

    for(t in 1:K){
        Rm[m == t] <- (1 - df$D) * (df$T == t)
        Lm[m == t] <- df$D * (df$T == t)
        Im[m == t] <- (df$T >= t)
        Jm[m == t] <- (df$T > t) * df$D + (df$T >= t) * (1 - df$D)
    }

    data <- data.frame(df[as.numeric(gl(n, K)), ], m = m, Im, Jm, Rm, Lm)

    return(data)

}

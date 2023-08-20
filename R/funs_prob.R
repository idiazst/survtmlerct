#' Compute TMLE of survival probability with efficiency gains
#'
#' @param data A long form data frame (e.g., the output of transformData*())
#' @param tau A time of interest
#'
#' @return Survival estimate for each treatment arm with influence function based standard error for the difference
#'
#' @examples
#' library(survival)
#' library(tidyverse)
#' library(dummies)
#'
#' data <- colon %>% dummy.data.frame(c('differ', 'extent')) %>%
#'   filter(rx != 'Obs') %>%
#'   mutate(A = rx == 'Lev+5FU', id = as.numeric(as.factor(id)),
#'          nanodes = is.na(nodes), nodes = ifelse(is.na(nodes), 0, nodes)) %>%
#'   select(-rx) %>%  group_by(id) %>% summarise_all(funs(min)) %>% select(-study) %>%
#'   rename(T = time, D = status)
#'
#' dlong <- transformData(data, 30)
#'
#' fitL <- glm(Lm ~ A * (m + sex + age + obstruct + perfor + adhere + nodes +
#'                         differ1 + differ2 + differ3 + differNA + extent1 + extent2 +
#'                         extent3 + extent4 + surg + node4 + etype),
#'             data = dlong, subset = Im == 1, family = binomial())
#' fitR <- glm(Rm ~ A * (as.factor(m) + sex + age + obstruct + perfor + adhere + nodes +
#'                         differ1 + differ2 + differ3 + differNA + extent1 + extent2 +
#'                         extent3 + extent4 + surg + node4 + etype),
#'             data = dlong, subset = Jm == 1, family = binomial())
#' fitA <- glm(A ~ sex + age + obstruct + perfor + adhere + nodes +
#'               differ1 + differ2 + differ3 + differNA + extent1 + extent2 +
#'               extent3 + extent4 + surg + node4 + etype,
#'             data = dlong, subset = m == 1, family = binomial())
#'
#' dlong <- mutate(dlong,
#'                 gR1 = bound01(predict(fitR, newdata = mutate(dlong, A = 1), type = 'response')),
#'                 gR0 = bound01(predict(fitR, newdata = mutate(dlong, A = 0), type = 'response')),
#'                 h1 = bound01(predict(fitL, newdata = mutate(dlong, A = 1), type = 'response')),
#'                 h0 = bound01(predict(fitL, newdata = mutate(dlong, A = 0), type = 'response')),
#'                 gA1 = bound01(predict(fitA, newdata = mutate(dlong, A = 1), type = 'response')))
#'
#' tau <- max(dlong$m)
#'
#' tmle_prob(dlong, tau)
tmle_prob <- function(data, tau){

    h1 <- data$h1
    h0 <- data$h0
    gR1 <- data$gR1
    gR0 <- data$gR0
    gA1 <- data$gA1[data$m == 1]
    id <- data$id
    A <- data$A

    n <- length(unique(id))
    m <- as.numeric(data$m)
    K <- max(m)

    gA0 <- 1 - gA1
    h  <- A*h1 + (1-A)*h0
    gR <- A*gR1 + (1-A)*gR0

    crit <- TRUE
    iter <- 1

    ind <- outer(m, 1:K, '<=')

    while(crit && iter <= 20){

        S1 <- tapply(1 - h1, id, cumprod, simplify = FALSE)
        S0 <- tapply(1 - h0, id, cumprod, simplify = FALSE)

        prodlag <- function(x) cumprod(c(1, x[-length(x)]))
        S1lag <- tapply(1 - h1, id, prodlag, simplify = FALSE)
        S0lag <- tapply(1 - h0, id, prodlag, simplify = FALSE)

        G1 <- tapply(1 - gR1, id, cumprod, simplify = FALSE)
        G0 <- tapply(1 - gR0, id, cumprod, simplify = FALSE)

        St1 <- do.call('rbind', S1[id])
        St0 <- do.call('rbind', S0[id])

        Sm1 <- unlist(S1)
        Sm0 <- unlist(S0)
        Sm1lag <- unlist(S1lag)
        Sm0lag <- unlist(S0lag)

        Gm1 <- unlist(G1)
        Gm0 <- unlist(G0)

        Z1 <- - (ind * St1)[, tau] / bound(Sm1 * gA1[id] * Gm1)
        Z0 <- - (ind * St0)[, tau] / bound(Sm0 * gA0[id] * Gm0)

        M <- tapply(Sm1 / bound(gA1[id]) * (m == tau), id, sum) +
            tapply(Sm0 / bound(gA0[id]) * (m == tau), id, sum)

        H1 <- - St1[, tau] / bound(Sm1lag * gA1[id] * Gm1)
        H0 <- - St0[, tau] / bound(Sm0lag * gA0[id] * Gm0)


        H <- A * H1 - (1-A) * H0

        eps   <- coef(glm2(Lm ~ 0 + offset(qlogis(h)) + I(A * Z1) + I((1-A) * Z0),
                           family = binomial(), subset = Im == 1, data = data))
        gamma <- coef(glm2(Rm ~ 0 + offset(qlogis(gR)) + H, family = binomial(),
                           subset = Jm == 1, data = data))
        nu    <- coef(glm2(A[m == 1] ~ 0 + offset(qlogis(gA1)) + M,
                           family = binomial()))

        h1old <- h1
        h0old <- h0
        gR1old <- gR1
        gR0old <- gR0
        gA1old <- gA1

        eps[is.na(eps)] <- 0
        gamma[is.na(gamma)] <- 0
        nu[is.na(nu)] <- 0

        h1 <- bound01(plogis(qlogis(h1) + eps[1] * Z1))
        h0 <- bound01(plogis(qlogis(h0) + eps[2] * Z0))
        gR1 <- bound01(plogis(qlogis(gR1) + gamma * H1))
        gR0 <- bound01(plogis(qlogis(gR0) - gamma * H0))
        gA1 <- bound01(plogis(qlogis(gA1) + nu * M))
        gA0 <- 1 - gA1
        h  <- A * h1  + (1 - A) * h0
        gR <- A * gR1 + (1 - A) * gR0
        iter <-  iter + 1

        crit <- any(abs(c(eps, gamma, nu)) > 1e-3/n^(0.6))

    }

    S1 <- tapply(1 - h1, id, cumprod, simplify = FALSE)
    S0 <- tapply(1 - h0, id, cumprod, simplify = FALSE)

    G1 <- tapply(1 - gR1, id, cumprod, simplify = FALSE)
    G0 <- tapply(1 - gR0, id, cumprod, simplify = FALSE)

    St1 <- do.call('rbind', S1[id])
    St0 <- do.call('rbind', S0[id])

    Sm1 <- unlist(S1)
    Sm0 <- unlist(S0)
    Gm1 <- unlist(G1)
    Gm0 <- unlist(G0)

    Z1 <- - (ind * St1)[, tau] / bound(Sm1 * gA1[id] * Gm1)
    Z0 <- - (ind * St0)[, tau] / bound(Sm0 * gA0[id] * Gm0)


    DT <- with(data, tapply(Im * (A * Z1 - (1 - A) * Z0) * (Lm - h), id, sum))
    DW1 <- with(data, St1[m == 1, tau])
    DW0 <- with(data, St0[m == 1, tau])
    theta1 <- mean(DW1)
    theta0 <- mean(DW0)
    theta <- theta1 - theta0
    D <- DT + DW1 - DW0 - theta
    sdn <- sqrt(var(D) / n)

    out <- list(prob = c(prob0 = theta0, prob1 = theta1), std.error.diff = sdn)
    class(out) <- 'adjusted_prob'
    return(out)

}

#' Compute AIPW of survival probability
#'
#' @param data A long form data frame (e.g., the output of transformData*())
#' @param tau A time of interest
#'
#' @return survival probability estimate for each treatment arm with influence function based standard error for the difference
#'
#' @examples
#' # See example of tmle() for the creation of dlong
#' aipw_prob(dlong, tau)
aipw_prob <- function(data, tau){

    h1 <- data$h1
    h0 <- data$h0
    gR1 <- data$gR1
    gR0 <- data$gR0
    gA1 <- data$gA1[data$m == 1]
    id <- data$id
    A <- data$A

    n <- length(unique(id))
    m <- as.numeric(data$m)
    K <- max(m)

    gA0 <- 1 - gA1
    h  <- A * h1 + (1 - A) * h0
    gR <- A * gR1 + (1 - A) * gR0

    ind <- outer(m, 1:K, '<=')

    S1 <- tapply(1 - h1, id, cumprod, simplify = FALSE)
    S0 <- tapply(1 - h0, id, cumprod, simplify = FALSE)

    G1 <- tapply(1 - gR1, id, cumprod, simplify = FALSE)
    G0 <- tapply(1 - gR0, id, cumprod, simplify = FALSE)

    St1 <- do.call('rbind', S1[id])
    St0 <- do.call('rbind', S0[id])

    Sm1 <- unlist(S1)
    Sm0 <- unlist(S0)
    Gm1 <- unlist(G1)
    Gm0 <- unlist(G0)

    Z1 <- -(ind * St1)[, tau] / bound(Sm1 * gA1[id] * Gm1)
    Z0 <- -(ind * St0)[, tau] / bound(Sm0 * gA0[id] * Gm0)

    DT1 <- with(data, tapply(Im * A * Z1 * (Lm - h), id, sum))
    DT0 <- with(data, tapply(Im * (1 - A) * Z0 * (Lm - h), id, sum))
    DW1 <- with(data, St1[as.numeric(m) == 1, tau])
    DW0 <- with(data, St0[as.numeric(m) == 1, tau])
    aipw <- c(mean(DT0 + DW0), mean(DT1 + DW1))

    D <- DT1 - DT0 + DW1 - DW0
    sdn <- sqrt(var(D) / n)
    
    out <- list(prob = c(prob0 = aipw[1], prob1 = aipw[2]), std.error.diff = sdn)
    class(out) <- 'adjusted_prob'
    return(out)
}

#' Compute unadjusted Kaplan-Meier estimators of survival probability
#'
#' @param data A long form data frame (e.g., the output of transformData*())
#' @param tau A time of interest
#'
#' @return survival probability estimate for each treatment arm
#'
#' @examples
#' # See example of tmle() for the creation of dlong
#' unadjusted(dlong, tau)
unadjusted_prob <- function(data, tau){

    fith  <- lm(Lm ~ as.factor(m) * as.factor(A), subset = Im == 1, data = data)
    fitgR <- lm(Rm ~ as.factor(m) * as.factor(A), subset = Jm == 1, data = data)
    fitgA <- lm(A ~ 1, subset = m == 1, data = data)

    h1 <- bound01(predict(fith, newdata = data.frame(m = data$m, A = 1),
                          type = 'response'))
    h0 <- bound01(predict(fith, newdata = data.frame(m = data$m, A = 0),
                          type = 'response'))
    gR1 <- bound01(predict(fitgR, newdata = data.frame(m = data$m, A = 1),
                           type = 'response'))
    gR0 <- bound01(predict(fitgR, newdata = data.frame(m = data$m, A = 0),
                           type = 'response'))
    gA1 <- bound01(predict(fitgA, type = 'response'))

    id <- data$id
    n <- length(unique(id))
    m <- as.numeric(data$m)
    K <- max(m)
    A <- data$A

    gA0 <- 1 - gA1
    h  <- A * h1 + (1 - A) * h0
    gR <- A * gR1 + (1 - A) * gR0
    ind <- outer(m, 1:K, '<=')
    S1 <- tapply(1 - h1, id, cumprod, simplify = FALSE)
    S0 <- tapply(1 - h0, id, cumprod, simplify = FALSE)
    G1 <- tapply(1 - gR1, id, cumprod, simplify = FALSE)
    G0 <- tapply(1 - gR0, id, cumprod, simplify = FALSE)
    St1 <- do.call('rbind', S1[id])
    St0 <- do.call('rbind', S0[id])
    Sm1 <- unlist(S1)
    Sm0 <- unlist(S0)
    Gm1 <- unlist(G1)
    Gm0 <- unlist(G0)

    DW1 <- with(data, St1[m == 1, tau])
    DW0 <- with(data, St0[m == 1, tau])
    km <- c(mean(DW0), mean(DW1))

    Z1 <- -(ind * St1)[, tau] / bound(Sm1 * gA1[id] * Gm1)
    Z0 <- -(ind * St0)[, tau] / bound(Sm0 * gA0[id] * Gm0)

    DT1 <- with(data, tapply(Im * A * Z1 * (Lm - h), id, sum))
    DT0 <- with(data, tapply(Im * (1 - A) * Z0 * (Lm - h), id, sum))
    D <- DT1 - DT0 + DW1 - DW0
    sekm <- sqrt(var(D) / n)

    out <- list(prob = c(prob0 = km[1], prob1 = km[2]), std.error.diff = sekm)
    class(out) <- 'adjusted_prob'
    return(out)

}

###########################################################
#### Stochastic epidemics models with random duration #####
################ Functions for simulations ################
###########################################################


### Number of contacts N(t)
# Values of N(t)
nk <- function(N.distr, N.param1){
    if (N.distr == "constant") {
        N.nk <- N.param1
    } else if (N.distr == "binomial" || N.distr == "trunc.poisson") {
        N.nk <- 0:N.param1
    } else if (N.distr == "neg.binomial") {
        N.nk <- 0:100       # OBS: It is fixed!
    }
    N.nk
}

# Density function of N(t)
pkt <- function(values, t, N.distr, N.param1, N.param2 = NULL){
    if (N.distr == "constant") {
        N.pkt <- 1
    } else if (N.distr == "binomial") {
        N.pkt <- dbinom(values, size = N.param1, prob = N.param2(t))
    } else if (N.distr == "trunc.poisson") {
        N.pkt <- dtpois(values, lambda = N.param2(t), b = N.param1)
    } else if (N.distr == "neg.binomial") {
        N.pkt <- dnbinom(values, size = N.param1, prob = N.param2)
    } 
    N.pkt
}



### Duration probabilities q_i
duration_prob <- function(q.distr, d, q.param = NULL){
    if (q.distr == "determinist") {
        qi <- NULL
    } else if (q.distr == "binomial") {
        qi <- 1 - pbinom(0:(d-1), size = d, prob = q.param)
    } else if (q.distr == "trunc.geom") {
        qi <- 1 - pgeom(0:(d-1), prob = q.param) / pgeom(d, prob = q.param)
    } else if (q.distr == "trunc.pois") {
        qi <- 1 - ptpois(0:(d-1), lambda = q.param, b = d)
    } else if (q.distr == "HKstudy") {
        qi <- c(0.857, 0.708, 0.56, 0.411, 0.31, 0.22, 0.071, 0.03, 0.006)
    }
    qi
}


### Transition probability
# Probability that a susceptible becomes infective
pty <- function(n, y, p, lambda, t, N.distr, N.param1, N.param2 = NULL){
    N.nk <- nk(N.distr, N.param1)
    N.pkt <- pkt(values = N.nk, t, N.distr, N.param1, N.param2)
    
    if (y == n-1) {
        prob <- 1-sum((1-p)^N.nk*N.pkt)
    } else {
        prob <- 1-sum((1-p*y*lambda/(n-1))^N.nk*N.pkt)
    }
    prob
}





### Model

## Notation
# Y     = Number of infectives
# Ynew  = Number of new infectives
# X     = Number of susceptibles

## Arguments for function Model()
# n         = Population size.
# Y0        = Number of infectives at time t0.
# p         = Probability of transmission.
# r         = Minimum duration of the disease.
# d         = "Stochastic" duration of the disease.
# q.distr   = Distribution of the duration probabilities q_i. It can take values: "determinist" or "binomial".
# q.param...= Parameters of the distribution of q_i:
#               "constant": q.param = NULL (d=0).
#               "binomial": d = size, q.param = prob.
# lambda    = Quarantine. It can take values in [0,1].
# Time      = Total time (days).
# N.distr   = Distribution of N(t). It can take values: "constant", "binomial" or "trunc.poisson".
# N.param...= Parameters of the distribution of N(t):
#               "constant": N.param1 = N.
#               "binomial": N.param1 = size, N.param2(t) = prob (depends on t).
#               "trunc.poisson": N.param1 = upper value (b of truncated poisson), N.param2(t) = lambda (depends on t).
#               "neg.binomial": N.param1 = size (dispersion parameter r, can be not integer), N.param2 = prob.

## Output of function Model(). List of:
# time      = The time the disease ends.
# Ytot      = Total number of infected.
# Y         = Vector of number of infectives at time t. 
# Ynew      = Vector of number of new infectives at time t.
# R0        = Value of R0. If Y0!=1, return NULL.

Model <- function(n, Y0, p, r, d = 0, q.distr, q.param = NULL, lambda, Time, N.distr, N.param1, N.param2 = NULL){
    ## Definition of internal variables
    Ynew <- Ynewt <- Y0
    Y <- Yt <- Y0
    X <- Xt <- n-Y0
    # Duration probabilities q_i
    if (q.distr == "determinist") {d <- 0}
    qi <- duration_prob(q.distr = q.distr, d = d, q.param = q.param)
    
    ## Model of the epidemic
    time <- 1
    while (Yt != 0 && time < Time) {
        # New infecteds
        Ynewt <- rbinom(1, Xt, pty(n = n, y = Yt, p = p, lambda = lambda, t = time-1, N.distr = N.distr, N.param1 = N.param1, N.param2 = N.param2))
        Ynew <- c(Ynew, Ynewt)
        
        # Infecteds
        Ynew_stillinf <- Ynew[max(1,time-r-d+2):(time+1)]
        q <- c(rep(1,r), qi)
        if (time < r+d-1) {q <- rev(q[1:(time+1)])} else {q <- rev(q)}
        Yt <- sum(mapply(rbinom, n = 1, size = Ynew_stillinf, prob = q))

        Y <- c(Y, Yt)
        
        # Susceptibles
        Xt <- n-sum(Ynew)
        X <- c(X, Xt)
        
        time <- time+1
    }
    
    # Calculation of R0
    if (Y0 == 1) {
        R <- r + rbinom(1, d, 0.5)
        R0 <- sum(Ynew[2:(R+1)] / Y[1:R])
    } else {
        R0 <- NULL
    }
    
    # Output
    infectives <- list(time = time-1, Ytot = sum(Ynew), Y = Y, Ynew = Ynew, R0 = R0)
    infectives
}


### Simulations

## Arguments for function simulations()
# n.trials  = Number of simulations.
# n         = Population size.
# Y0        = Number of infectives at time t0.
# p         = Probability of transmission.
# r         = Duration of the disease.
# d         = "Stochastic" duration of the disease.
# q.distr   = Distribution of the duration probabilities q_i. It can take values: "determinist" or "binomial".
# q.param...= Parameters of the distribution of q_i:
#               "constant": q.param = NULL (d=0).
#               "binomial": d = size, q.param = prob.
# lambda    = Quarantine. It can take values in [0,1].
# Time      = Total time (days).
# N.distr   = Distribution of N(t). It can take values: "constant", "binomial" or "trunc.poisson".
# N.param...= Parameters of the distribution of N(t):
#               "constant": N.param1 = N.
#               "binomial": N.param1 = size, N.param2(t) = prob (depends on t).
#               "trunc.poisson": N.param1 = upper value (b of truncated poisson), N.param2(t) = lambda (depends on t).
#               "neg.binomial": N.param1 = size (dispersion parameter r, can be not integer), N.param2 = prob.

## Output of function simulations(). List of:
# time      = Vector of the times the disease ends.
# Ytot      = Vector of the total number of infected.
# Y         = Matrix of number of infectives at time t. Each column is a simulation. 
# Ynew      = Matrix of number of new infectives at time t. Each column is a simulation. 
# R0        = Vector of R0. If Y0!=1, return NULL.

simulations <- function(n.trials, n, Y0, p, r, d = 0, q.distr, q.param = NULL, lambda, Time, N.distr, N.param1, N.param2 = NULL) {
    # Description of the simulations
    cat("This function runs", n.trials, "simulations of the spread of a disease with the following properties: \n",
        "Population size:", n,
        "\n Infected at time 0:", Y0,
        "\n Probability of transmission:", p,
        "\n Duration of the disease:", q.distr, "with minimum duration", r, "and stochastic duration", d, "with parameters", q.param,
        "\n Immunity parameter (lambda):", lambda,
        "\n Total time considered:", Time,
        "\n Distribution of the contacts:", N.distr, "with parameters", N.param1, "and", substitute(N.param2), "\n")

    # Simulations
    simul <- replicate(n = n.trials, Model(n, Y0, p, r, d, q.distr, q.param, lambda, Time, N.distr, N.param1, N.param2), simplify = FALSE)
    
    # Outputs
    time <- sapply(simul, "[[", 1)
    Ytot <- sapply(simul, "[[", 2)
    
    Y <- lapply(simul, "[[", 3)
    Y <- sapply(Y, "[", 1:max(time))
    Y[is.na(Y)] <- 0
    
    Ynew <- lapply(simul, "[[", 4)
    Ynew <- sapply(Ynew, "[", 1:max(time))
    Ynew[is.na(Ynew)] <- 0
    
    R0 <- sapply(simul, "[[", 5)
    
    # Outputs list
    simul.list <- list(time = time, Ytot = Ytot, Y = Y, Ynew = Ynew, R0 = R0)
    simul.list
}


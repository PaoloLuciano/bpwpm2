# Math Utils for the bpwpm2 package
#-------------------------------------------------------------------------------

#' Piece wise polinomial expansion for X (PWP)
#'
#' Calculates and returns a new designer matrix, representing the correspondieng PWP
#' expansion for all d variables. Combines all the parameters on a relatively
#' fast computation of basis expansion for X described on the thesis and on
#' (Denison, Mallik and Smith, 1998).
#'
#' @inheritParams bpwpm2_gibbs
#' @param n
#' @param d
#' @param tau
#'
#' @return A desginer matrix containing the PWP expansion described in (...)
#' @export
#'
calculate_Psi <- function(X, M, J, K, n ,d, tau){

    # Psi is a list containing the basis transformations matrices for each dimension j.
    # For now, this basis expansion is done following the formula on the thesis, optimized as much as posible

    # A diagram of this expansion can be found on the thesis
    Psi <- rep(1, times = n)

    for(j in seq(1,d)){

        # Creating the first basis polinomial
        Psi_partial <- sapply(X = seq(1, M - 1), FUN = function(x , w){w ^ x}, w = X[ , j])

        # Piecewise part
        for(k in seq(1,J-1)){
            Psi_partial <- cbind(Psi_partial, sapply(X = seq(K, M-1),
                                                             FUN = function(x,w){w^x},
                                                             w = pmax(0, X[ ,j] - tau[k,j])))
        }
        Psi <- cbind(Psi,Psi_partial)
    }

    return(Psi)
}

#-------------------------------------------------------------------------------

#' F matrix calculation
#'
#' Function for calculating the F matrix, described on the thesis.
#' This is the transformed input matrix that depends on the piecewise polinomial
#' expansion Psi and a set of weights w.
#'
#' @param Psi Piecewise Polinomail expansion for an input matrix X previosly
#'  calculated by \code{\link{calculate_Psi}}
#' @param beta Posterior estimate of the beta parameters previously calculated
#' with \code{\link{posterior_params}}
#' @param d Number of dimentions, this parameter helps to improve efficiency
#'
#' @return F matrix
#' @export
#'
calculate_F <- function(Psi, beta, d){

    dim_Psi <- dim(Psi)
    n <- dim_Psi[1]
    lambda <- dim_Psi[2] # Total number of parameters
    N <- (lambda - 1)/d

    # Index to split the matrixes and vectors
    index <- seq(2, length.out = N)

    # Calculating F matrix
    # F_0
    mat_F <- rep(beta[1], n)

    for(j in seq(1,d)){
        mat_F <- cbind(mat_F, crossprod(t(Psi[ , index]),beta[index]))
        index <- index + N
    }

    return(mat_F)
}

#-------------------------------------------------------------------------------
#' Ergodic Mean
#'
#' Calculates the Ergodic Mean for an MCMC_Chain Matrix producede as output from
#' \code{\link{bpwpm_chain}}. It is used by the \code{\link{plot_ergodic_mean}}
#'  function, but left available for the user.
#' @param mcmc_chain
#'
#' @return The Ergodic Mean Matrix
#' @export
#'
#' @examples MA_betas <- (betas)
ergodic_mean <- function(mcmc_chain){
    return(apply(mcmc_chain,2,cumsum)/seq(1:dim(mcmc_chain)[1]))
}

#-------------------------------------------------------------------------------

#' Log Loss
#'
#' An implementation of the Log-Loss function for the binomial case
#'
#' @inheritParams bpwpm_gibbs
#' @param p Vector of fitted probabilities for each Y.
#' @param eps Machine error to hanlde limit cases on the logarithm function
#'
#' @return The value of the Log Loss function (numeric). The smaller, the better
#' @export
#'
#' @examples log_loss(true_values, fitted probabilities)
#' log_loss(true_values, fitted probabilities)
#' log_loss(true_values, fitted probabilities, 1e-30, FALSE)
log_loss <- function(Y, p, eps = 1e-10, verb = TRUE){

    if(class(Y) == "factor"){
        Y <- as.integer(Y) -1
    }

    p_corrected = pmin(pmax(p, eps), 1-eps)
    ll <- - sum (Y * log(p_corrected) + (1 - Y) * log(1 - p_corrected))/length(Y)

    if(verb){
        cat("\nLog_Loss: ",ll, sep ="")
    }

    return(- sum (Y * log(p_corrected) + (1 - Y) * log(1 - p_corrected))/length(Y))
}

#-------------------------------------------------------------------------------

#' Calculate the Accuracy of the model
#'
#' Given a set of true values and their corresponding fitted probabilities,
#' the function calculates the accuracy of the model defined by:
#'  \eqn{1- #wrong prediction/# of observations}
#'
#' @inheritParams  log_loss
#'
#' @return The accuracy of the model, given the fitted probabilities and new data
#' @export
#'
#' @examples (new_Y, fitted_probs_for_data)
accuracy <- function(new_Y, p, verb = FALSE){

    if(class(new_Y) == "factor"){
        new_Y <- as.integer(Y) - 1
    }

    n <- length(new_Y)

    est_Y <- rep(0, times = n)
    est_Y[p > 0.5] <- 1

    wrongs <- sum(abs(new_Y - est_Y))

    if(verb) cat(wrongs, " incorrect categorizations\n", sep = "")

    return(1 - wrongs/n)
}

#-------------------------------------------------------------------------------
#' Contingency Table for the prediciton of a bpwpm
#'
#' @param new_Y Response variable to test the model for
#' @param p Vector of fitted probabilities
#'
#' @return The contingency table
#' @export
#'
#' @examples contingency_table(train_y, est_p)
contingency_table <- function(new_Y, p){

    if(class(new_Y) == "factor"){
        new_Y <- as.integer(Y) - 1
    }

    n <- length(new_Y)
    est_Y <- rep(0, times = n)
    est_Y[p > 0.5] <- 1

    a <- sum((new_Y + est_Y) == 0)
    b <- sum((est_Y - new_Y == 1))
    c <- sum((new_Y - est_Y == 1))
    d <- sum((new_Y + est_Y) == 2)

    ct <- base::data.frame(rbind(
        cbind(a, b, sum(-(new_Y - 1))),
        cbind(c, d, sum(new_Y)),
        cbind(sum(-(est_Y - 1)), sum(est_Y), n)))

    colnames(ct) <- c("Est. Y = 0", "Est Y = 1", "Real Y - Totals")
    row.names(ct) <- c("Real Y = 0", "Real Y = 1", "Est. Y - Totals")

    return(ct)
}

#-------------------------------------------------------------------------------

#' Calculate prediction eta
#'
#' @param F_mat The PWP transformed input space, calculated by
#' \code{\link{calculate_F}}
#'
#' @return A numeric vector representing the eta of \code{R^d} into
#' \code{R} given all of the parametres
#' @export
#'
calculate_eta <- function(F_mat){
    # Since the model is additive and it doesn't depend on other parameters (thank god)
    # eta is just the sum of each F row
    return(apply(F_mat,1,sum))
}

#-------------------------------------------------------------------------------

#' Calculates trained model for new Data
#'
#' Given a set of parameters inherited by the \code{\link{bpwpm_gibbs}} training
#' procedure, we calculate the value of the eta function for new data.
#' This function is both used on the predict and plot_3D functions.
#' @param new_X A new set of data for which to calculate the f(x) eta
#' function
#' @param bpwpm_params A list of bpwpm parameters created by the function
#' \code{\link{posterior_params}}
#'
#' @return A vector of the eta vector for a given set of points
#' @export
#'
#' @examples (test_data, mean_params)
model_eta <- function(new_X, bpwpm_params){

    if(class(bpwpm_params) != "bpwpm_params"){
        error("Invalid class, object should be of class bpwpm_params")
        geterrmessage()
    }

    Psi <- calculate_Psi(X = new_X,
                         M = bpwpm_params$M,
                         J = bpwpm_params$J,
                         K = bpwpm_params$K,
                         n = dim(new_X)[1],
                         d = bpwpm_params$d,
                         tau = bpwpm_params$tau)

    F_mat <- calculate_F(Psi = Psi,beta = bpwpm_params$beta, d = bpwpm_params$d)

    return(calculate_eta(F_mat = F_mat))
}

#-------------------------------------------------------------------------------

#' Calculate Probabilities of Binary Outcome
#'
#' Given a model, we can calculate the corresponding fitted probabilites of the
#' random response variable Y. Whereas this is new data or the one used to train
#' the model. Since the model is a probit GLM at this point, we only need to
#' calculate the eta and then plug them on the inverse of the normal
#' gaussian acomulation function
#'
#' @inheritParams model_eta
#'
#' @return A vector of fitted probabilities for the response variable Y
#' @export
#'
posterior_probs <- function(new_X, bpwpm_params){
    if(class(bpwpm_params) != "bpwpm_params"){
        error("Invalid bpwpm parameters")
        geterrmessage()
    }

    eta <- model_eta(new_X,
                          bpwpm_params)

    return(pnorm(eta))
}









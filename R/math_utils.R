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

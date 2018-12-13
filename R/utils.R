# Series of utility functions to evaluate the MCMC chain for the bpwpm2 model
#-------------------------------------------------------------------------------

#' Thin MCMC Chain
#'
#' Thins and eliminates the burn in period of the Gibbs Chain.
#'
#' @param mcmc_chain An MCMC Chain matrix. (draws * number of params)
#' @inheritParams posterior_params
#'
#' @return The thinned down version of the MCMC chain
#' @export
#'
#' @examples (beta, 2, 2000)
thin_chain <- function(mcmc_chain, burn_in, thin){

    draws <- dim(mcmc_chain)[1]

    if(burn_in > draws){stop("Burn-in parameter too large")
        geterrmessage()}

    # From the burn_in parameter, return every (thin+1) index
    return(mcmc_chain[burn_in + seq(1, floor((draws - burn_in)/(thin + 1))), ])
}

#-------------------------------------------------------------------------------

#' Thinning of a BPWPM
#'
#' Thins the beta chain of an bpwpm2 and returns an object of the same
#' type
#'
#' @param bpwpm A bpwpm object created by \code{\link{bpwpm_gibbs}}
#' @inheritParams posterior_params
#'
#' @return An object of the same kind of the bpwpm but thinned down
#' @export
#'
thin_bpwpm <- function(bpwpm, burn_in, thin){

    if(class(bpwpm) != 'bpwpm'){
        error("Object not of class bpwpm")
        geterrmessage()
    }

    bpwpm_model_copy <- bpwpm
    bpwpm_model_copy$beta <- thin_chain(bpwpm$beta, burn_in = burn_in, thin = thin)

    return(bpwpm_model_copy)
}

#-------------------------------------------------------------------------------

#' Puntual Posterior Estimation of bpwpm2 beta parameters
#'
#' Given a model output by the function \code{\link{bpwpm_gibbs}}, it thins
#' down the chain and makes puntual estimation for the parameters, given a type
#' of estimation. This parameter object can later be used to calculate other
#' results.
#'
#' @param bpwpm an object of the class bpwpm
#' @param thin A thinning parameter for the MCMC Chain
#' @param burn_in A burn in parameter for the MCMC Chain
#' @param type The type of punctual estimation for the parameters. Options
#' include: mean, mode or median.
#'
#' @return An object of the class bpwpm_params that contains:
#' \describe{
#'   \item{\eqn{\beta}}{The posterior estimation for \eqn{\beta}}
#'   \item{\eqn{\tau}}{The corresponding nodes}
#'   \item{F}{The final F matrix}
#'   \item{M}{M parameter}
#'   \item{J}{J parameter}
#'   \item{K}{K parameter}
#'   \item{d}{number of dimentions}
#'   \item{indep_terms}{Logical}
#' }
#' @export
#'
#' @examples (model, 2, 2000, 'mean') (model, 0, 0, 'median')
posterior_params <- function(bpwpm, burn_in, thin, type = 'mean'){

    if(class(bpwpm) != 'bpwpm'){
        error("Object not of class bpwpm")
        geterrmessage()
    }

    if(type == 'mean'){func <- mean}
    else if(type == 'median'){func <- median}
    else{error("Incorrect type of parameter estimation")
        geterrmessage()}

    # Fist the model is thinned down
    thined_model <- thin_bpwpm(bpwpm, burn_in = burn_in, thin = thin)

    # Estimation for beta
    estimated_beta <- sapply(thined_model$beta, func)

    # Estimation of F matrix
    estimated_F <- calculate_F(bpwpm$Psi, estimated_beta, d = bpwpm$d)
    rownames(estimated_F) <- NULL
    colnames(estimated_F) <- paste("f_", seq(0, bpwpm$d), sep = "")


    params <- list(beta = estimated_beta,
                   tau = bpwpm$tau,
                   estimated_F = estimated_F,
                   M = bpwpm$M,
                   J = bpwpm$J,
                   K = bpwpm$K,
                   d = bpwpm$d,
                   indep_terms = bpwpm$indep_terms)
    class(params) <- 'bpwpm_params'

    return(params)
}

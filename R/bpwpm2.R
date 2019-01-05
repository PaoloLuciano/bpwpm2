#' Bayesian Piece-Wise Polynomial Model 2 Algorithm
#'
#' Second take at the bpwpm model, based solely on the Albert + Chibb Paper
#'
#' @param Y Response vector of n binary observatios (integers 0,1 - vector of
#'   size n) Can be encoded as a factor a numeric vector.
#' @param X Design matrix of n observations and d covariables (numeric - n*d)
#' @param M M minus 1 is the degree of the polinomial (integer - M > 0)
#' @param J Number of intervals in each dimention (integer - J > 1)
#' @param K Order of continuity in the derivatives (integrer - 0 < K < M)
#' @param draws Númber of samples to draw from the Gibbs Sampler (integer - draw
#'   > 0)
#' @param tau the initial position of nodes selected by the user. although·
#' arbitraty they need to match the dimentions. (numeric - (J-1)*d)
#' @param beta_init Inital value for the Gibbs Sampler Chain (numeric - vector of
#'   size 1*(1 + Nd))
#' @param mu_beta_0 Prior Mean of Beta (numeric - matrix of size N*d)
#' @param sigma_beta_0_inv sigma_w_0_inv:= Prior Inverse if the Variance-Covariance
#'   Matrices of w (list - d elements, each element is a numeric matrix of size
#'   N*N)
#' @param precision_w If using the default sigmas for w, a diagonal matrix will
#'   be used. Precision controls its magnitude (numeric - precision > 0)
#' @param eps Numerical threshold
#' @param verb short for verbose, if TRUE, prints aditional information (logical)
#' @param debug_verb If TRUE, print even more info to help with debug_verbging (logical)
#'
#' @return An object of the class "bpwpm" containing at the following
#' components:
#' \describe{
#' \item{beta: }{A data frame containing the Gibbs sampler simulation for beta}
#' \item{Psi: }{The PWP Expansion for input matrix X and nodes selected on
#' percentiles}
#' \item{tau: }{Nodes used for training}
#' \item{M: }{Initial parameters}
#' \item{J: }{Initial parameters}
#' \item{K: }{Initial parameters}
#' \item{d: }{Number of dimentions}
#' \item{indep_terms: }{Logical. If independent terms are keept}
#' \item{info}{A string that prints the basic information of the mode. Used for
#' the summary function.}
#' }
#' @export
#'
#' @examples See the main document of thesis for a couple of full examples
#' with its corresponding analysis.
bpwpm2_gibbs <- function(Y, X, M, J, K,
                        draws = 10^3,
                        tau = NULL,
                        beta_init = NULL,
                        mu_beta_0 = NULL, sigma_beta_0_inv = NULL, precision_beta = 1,
                        eps = 1e-15,
                        verb = FALSE, debug_verb = FALSE){

    # 0. Stage Setting----------------------------------------------------------

    # 0.1 initializing parameters
    dim_X <- dim(X)
    n <- dim_X[1] # Number of observations
    d <- dim_X[2] # Number of covariables

    N <- M * J - K * (J - 1) - 1 # Number of basis funcitons for each d
    lambda <- 1 + d*N # Total number of parameters

    if(is.null(beta_init)){
        beta <- rep(1, times = lambda)
        temp_string <- paste("Utilizing standard initializing parameters")
    }else{
        beta <- beta_init
        temp_string <- paste("Utilizing custom parameters")
    }

    if(is.null(mu_beta_0)){
        mu_beta_0 <- rep(0, times = lambda)
    }

    if(is.null(sigma_beta_0_inv)){
        sigma_beta_0_inv <- precision_beta * diag(lambda)
    }

    # O.2 Dimension and parameters check
    if(!all(identical(dim_X[1],length(Y)), K < M, M > 1, J > 1, K > 0, n > lambda)){
        stop("Error in dimensionalities or parameters")
        geterrmessage()
    }

    rm(dim_X)

    # 0.3 Sanitizing Inputs
    if(class(Y) == "factor"){
        Y <- as.integer(Y) - 1
    }

    # 0.4 Basic Info print
    info <- paste("\n\tBPWPM MODEL\n
                  \tDimensions and pararamters check\n\t",
                  "Algorithm based on d = ", d, " covariables", "\n\t",
                  "Number of nodes J - 1 = ", J - 1, "\n\t",
                  "Order of polinomial M - 1 = ", M - 1, "\n\t",
                  "Order of continuity on derivatives K = ", K, "\n\t",
                  "Number of basis functions N = ", N, "\n\t",
                  "Number of observations n = ", n, "\n\t",
                  "Total number of params = ", lambda, "\n\t",
                  temp_string, "\n\n", sep = "")
    cat(info)

    # 1. Initializing Gibbs Sampler---------------------------------------------

    cat("Initializing Gibbs Sampler\n")

    # 1.1 Node Initialization.
    # Setting Nodes on the quantiles of X. (numeric, (J-1)*d)
    if(is.null(tau)){
        tau <- sapply(data.frame(X), quantile, probs = seq(0,1, by = 1/J))
        tau <- matrix(tau[-c(1,J+1), ], nrow = J - 1, ncol = d)
    }else if(!(dim(tau)[1] == (J-1) && dim(tau)[2] == d)){
        error("Dimentions of the given tau does not match up. The standard tau matrix is recomended")
        geterrmessage()
    }

    if(verb){
        # Nodes
        cat("\tNodes Locations\n")
        print(tau, digits = 3)
        cat("\n")
    }

    # 1.2 Calculates the designer matrix based on the specified parameters and splits it for other purposes
    Psi <- calculate_Psi(X, M, J, K, n, d, tau)

    #---------------------------- Matrix split ------------------------------------------------------------

    # 1.3 Final variable initialization
    sim_beta <- list()

    # For speeding up calculations
    Psi_cross <- crossprod(Psi,Psi)
    sigma_beta <- solve(sigma_beta_0_inv + Psi_cross) # Invariant throught the calculations
    z <- rep(0, times = n)
    counter_print <- 10^(ceiling(log10(draws)) - 1)

    # 2. Gibbs Sampler, z -> beta -> z ... ------------------------------
    for(k in seq(1,draws)){

        # Print number of iterations
        if(debug_verb){
            cat("\nIter: ", k, sep ="")
        }else if((k %% counter_print) == 0){
            cat("\nIter: ", k, sep = "")
        }

        # Linear predictor
        eta <- crossprod(t(Psi),beta)

        # 2.1. Z - Sampling from the truncated normal distribution for Z, taking into account the numerical erros
        temp_probs <- pnorm(0,eta,1)
        temp_probs = pmin(pmax(temp_probs, eps), 1 - eps)

        z[Y == 0] <- qnorm(runif(n = sum(1-Y), min = 0,
                                 max = temp_probs[Y == 0]), eta[Y == 0],1)
        z[Y == 1] <- qnorm(runif(n = sum(Y),
                                 min = temp_probs[Y == 1],max = 1), eta[Y == 1], 1)

        # 2.2. BETA - Sampling from the final distribution for beta
        mu_beta <- crossprod(t(sigma_beta),(crossprod(t(sigma_beta_0_inv),mu_beta_0) +
                                        crossprod(Psi,z)))
        beta <- c(mvtnorm::rmvnorm(1,mu_beta,sigma_beta)) # Simulating from the resulting distribution

        # Making beta simulation matrix
        if(k == 1){
            sim_beta <- beta # First iteration
        }else{
            sim_beta <- rbind(sim_beta,beta) # Apending Beta
        }

        if(verb){
            cat("\n\tBeta:\t", format(beta,digits = 2, width = 10), sep = "")
        }
    }

    # 3. Creating Output--------------------------------------------------------

    # 3.1 Naming Beta
    rownames(sim_beta) <- NULL
    colnames(sim_beta) <- paste("beta_", seq(0, lambda - 1), sep = "")

    model <- list(beta = data.frame(sim_beta),
                  Psi = Psi,
                  tau = tau,
                  M = M,
                  J = J,
                  K = K,
                  d = d,
                  info = info)
    class(model) <- 'bpwpm'

    beepr::beep(1) # Notifies when the simulation has finished
    return(model)
}

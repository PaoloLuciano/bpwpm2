# Plot Functionality for package bpwpm
#-------------------------------------------------------------------------------

#' Plot MCMC Chains
#'
#' Plots the last n draws of an MCMC chain
#'
#' @inheritParams plot.bpwpm
#' @param mcmc_chain An MCMC Chain matrix. (draws * number of params)
#' @param n Number of draws to plot
#' @param title Title for the plot
#'
#' @return A ggplot2 lines plot
#' @export
#'
#' @examples plot_chains(betas), plot_chains(w_j, 1000)
plot_chains <- function(mcmc_chain, n = 100, title = ""){

    dim_mcmc <- dim(mcmc_chain)
    n <- min(n, dim_mcmc[1])

    mcmc_temp <- tidyr::gather(mcmc_chain[seq(dim_mcmc[1] - n + 1,dim_mcmc[1]),],
                               key = Parameters)

    ggplot2::ggplot(mcmc_temp, aes(x = rep(seq(1,n),dim_mcmc[2]),
                                   y = value, group = Parameters,
                                   colour = Parameters)) +
        geom_line() + xlab("Index") + ylab("Value") + ggtitle(title)
}

#-------------------------------------------------------------------------------
#' Plot Ergodic Mean
#'
#' Plots the Ergodic Mean of an object of class \code{bpwpm}
#'
#' @inheritParams plot.bpwpm
#' @inheritParams thin_chain
#'
#' @return A series of plots for the ergodic mean of the parameters
#' @export
#'
#' @examples (model1, 0, 0)
plot_ergodic_mean <- function(object, thin = 0, burn_in = 0, ...){

    # if(!('bpwpm' %in% class(object))){
    #     error("Object not of the class bpwpm")
    #     geterrmessage()
    # }

    # thin_object <- thin_bpwpm(object, thin, burn_in)
    # Plots the whole erg
    # n <- dim(thin_object$betas)[1]
    # em_temp <- ergodic_mean(thin_object$betas)

    em_temp <- ergodic_mean(object)
    p <- plot_chains(data.frame(em_temp), title = "Betas - Ergodic Mean", ...)
    print(p)
}

#-------------------------------------------------------------------------------

#' Scatter Plot of 2D data
#'
#' Scatter plot to visualize the data and it's corresponding groups.
#'
#' @inheritParams plot_2D
#'
#' @return A ggplot2 scatter plot
#' @export
#'
#' @examples (Y = rbinom(100, 1, 4), X = cbind(rnorm(100), rnorm(100)))
plot_2D_data <- function(Y,X){

    if(dim(X)[2] != 2){
        error("Only a 2D plot can be made. X matrix has diferent dimensions")
        geterrmessage()
    }

    if(class(Y) == "numeric"){
        Y <- factor(Y)
    }

    if(class(X) == "matrix"){
        X <- data.frame(X)
    }

    ggplot2::ggplot(data = X, aes(x = X[, 1], y = X[ ,2], col = Y)) +
        geom_point() + xlab("X_1") + ylab("X_2")
}

#-------------------------------------------------------------------------------

#' Plot 2D projection of the Model
#'
#' 2D projection of both the inputs and the posterior regions of classification.
#' Usefull to evaluate the corresponding binary outcomes.
#' Instead of plotting the corresponding conotur of the 3D function plotted by
#' \code{\link{plot_3D_proj}}. The projection function is mapped to it's
#' corresponding binary output and plotted behind the regular data.
#'
#' @inheritParams plot_2D
#'
#' @return A ggplot2 scatter plot
#' @export
#'
plot_2D_proj <- function(Y, X, bpwpm_params, n = 15, alpha = 0.6){

    if(dim(X)[2] != 2){
        error("Only a 2D plot can be made. X matrix has diferent dimensions")
        geterrmessage()
    }

    if(class(Y) == "numeric"){
        Y <- factor(Y)
    }

    if(class(X) == "matrix"){
        X <- data.frame(X)
    }

    mins <- apply(X, 2, min)
    maxs <- apply(X, 2, max)

    linspace <- expand.grid(X1 = seq(mins[1] - 0.2, maxs[1] + 0.2, by = 1/n),
                            X2 = seq(mins[2] - 0.2, maxs[2] + 0.2, by = 1/n))

    linspace$Y <- model_eta(new_X = linspace,
                            bpwpm_params = bpwpm_params)

    linspace$Y<-  as.factor(as.integer(linspace$Y >= 0))
    linspace$a <- rep(alpha, times = dim(linspace)[1])

    data <- data.frame(cbind(X,Y), a = rep(1, times = dim(X)[1]))
    colnames(data) <- c("X1","X2","Y","a")

    data <- data.frame(rbind(data,linspace))

    ggplot2::ggplot(data = data) +
        geom_point(aes(x = X1, y = X2, col = Y, alpha = a), show.legend = FALSE) +
        xlab("X_1") + ylab("X_2")
}

#-------------------------------------------------------------------------------

#' Plots each dimention f(x)
#'
#' With the posterior \code{w} parameters calculated from the Gibbs run, and
#' \code{\link{posterior_params}}, a Final F matrix can be calculated. and
#' hence, ploted against every Input X to see how does the PWP expansion looks
#' like for the specified set of parameters.
#'
#' @param Y A vector of binary response. Can be encoded as either a factor
#' vector or as a numeric one.
#' @param X A data frame or matrix containing the original Inputs for the model.
#' @param F_mat The F matrix calculated via \code{\link{calculate_F}} or
#' alternativly you can pass it the parameters calculated by function
#' \code{\link{posterior_params}}
#'
#' @return d plots for each dimention created using ggplot2
#' @export
#'
plot_each_F <- function(Y, X, F_mat){

    if(class(F_mat) == "bpwpm_params"){

        if(dim(X)[1] == dim(F_mat$estimated_F)[1]){
            cat("Old F is being used")
            F_mat <- F_mat$estimated_F
        }
        else{
            cat("Calculating new F")
            M <- F_mat$M
            J <- F_mat$J
            K <- F_mat$K
            d <- F_mat$d
            tau <- F_mat$tau
            Phi <- calculate_Phi(X,M,J,K,d,tau, indep_terms = F_mat$indep_terms)
            F_mat <- calculate_F(Phi, F_mat$w, d)
        }
    }

    d <- dim(X)[2]

    if(class(Y) == "numeric" | class(Y) == "integer"){
        Y <- as.factor(Y)
    }

    for(i in seq(1:d)){
        p <- ggplot2::qplot(x = X[,i],  y = F_mat[,i+1],
                            color = Y) +
            xlab(paste("X_",i, sep = "")) +
            ylab(paste("F_",i,"(X_",i,")",sep = ""))
        print(p)
        if(i != (d+1)){
            readline(prompt="Press [enter] to view next plot")
        }
    }
}


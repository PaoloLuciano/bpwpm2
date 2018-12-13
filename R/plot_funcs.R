# Plot Functionality for package bpwpm
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#' Plot MCMC Chains
#'
#' Plots the last n draws of an MCMC chain
#'
#' @inheritParams plot.bpwpm
#' @param mcmc_chain An MCMC Chain matrix. (draws * number of params)
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

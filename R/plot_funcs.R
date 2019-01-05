# Plot Functionality for package bpwpm
#-------------------------------------------------------------------------------
#' Generic bpwpm object plotting
#'
#' Once a bpwpm has been run using the function \code{\link{bpwpm2_gibbs}}, the
#' chains can be plotted and hence evaluated. This generic function creates 3
#' set of graphs for the parameters \eqn{\beta} of the model:
#' - Ergodic means
#' - Trace of the parameters (in groups of size m)
#' - Histograms of the parameters (in groups of size m)
#'
#'
#' @param object of the class bpwpm
#' @param n number of draws to plot (the function plots from the end, ie:
#' n = 100 will plot the last 100 draws of the chain)
#' @param ... additional parameters to be passed to the functions
#'
#' @return a series of line plots and histograms
#' @export
#'
#' @examples (model1, 1000), (model2, 2000)
plot.bpwpm <- function(object, n = 1000, m = 5, ...){

    if(!('bpwpm' %in% class(object))){
        error("Object not of the class bpwpm")
        geterrmessage()
    }

    # Betas
    p <- plot_ergodic_mean(object$beta, n, title = paste("Betas: M = ", object$M,
                                                    ", J = ", object$J,
                                                    ", K = ", object$K, sep = ""))
    print(p)

    lambda <- dim(object$beta)[2] # Total number of parameters
    for(i in 1:ceiling(lambda/m)){
        readline(prompt = "Press [enter] to view next plot")

        index <- 1:m + m*(i-1)
        index <- index[index <= lambda]
        p <- plot_chains(object$beta[ ,index], n, title = paste("Betas: ",
                                                              min(index) - 1,
                                                              " to ",
                                                              max(index) - 1,
                                                              sep = ""))
        print(p)

        readline(prompt = "Press [enter] to view next plot")
        p <- plot_hist(object$beta[ ,index], n, title = paste("Betas: ",
                                                                min(index) - 1,
                                                                " to ",
                                                                max(index) - 1,
                                                                sep = ""))
        print(p)
    }
}

#-------------------------------------------------------------------------------
#' Plot Ergodic Mean
#'
#' Plots the Ergodic Mean of either a matrix of parameters or a bpwpm2 object.
#'
#' @inheritParams plot.bpwpm
#' @inheritParams thin_chain
#'
#' @return A series of plots for the ergodic mean of the parameters
#' @export
#'
#' @examples (model1, 0, 0)
plot_ergodic_mean <- function(object, n = 1000, thin = 0, burn_in = 0, ...){

    if(('bpwpm' %in% class(object))){
        # thin_object <- thin_bpwpm(object, thin, burn_in)
        # Plots the whole erg
        # n <- dim(thin_object$betas)[1]
        # em_temp <- ergodic_mean(thin_object$betas)

    }else{

    em_temp <- ergodic_mean(object)
    p <- plot_chains(data.frame(em_temp),n = n, ...)
    return(p)

    }
}

#'------------------------------------------------------------------------------
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
plot_chains <- function(mcmc_chain, n = 100, title = "", ...){

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
#' Plot MCMC Chains histograms
#'
#' Plot the parameters to test for convergence
#'
#' @inheritParams plot.bpwpm
#' @inheritParams plot_chains
#'
#' @return A histogram for the n draws and parameters of the chain
#' @export
#'
plot_hist <- function(mcmc_chain, number = 100, title = "", ...){

    dim_mcmc <- dim(mcmc_chain)
    n <- min(number, dim_mcmc[1])

    beta_temp <- tidyr::gather(mcmc_chain[seq(dim_mcmc[1] - n + 1,dim_mcmc[1]),],
                               key = Parameters)

    ggplot2::ggplot(beta_temp, aes(x = value, fill = Parameters)) +
        geom_histogram(..., position = "dodge") + xlab("Value") +
        ggtitle(title)

}
#-------------------------------------------------------------------------------

#' Generic function for plotting bpwpm_predictions objects
#'
#' Once a model has been run and evaluated, a prediction can be made using the
#' function \code{\link{predict.bpwpm}}. The input matrix \code{X} and
#' output \code{Y} are saved and can be plotted against the final PWP expansion
#' for the model.
#'
#' @param object of the class bpwpm_prediction
#' @param ... other arguments
#'
#' @return a series of plots from ggplot2
#' @export
#'
#' @examples (train_set_prediciton),
#'  (test_set_prediciton)
plot.bpwpm_prediction <- function(object, ...){

    if(!('bpwpm_prediction' %in% class(object))){
        error("Object not of the class bpwpm_prediction")
        geterrmessage()
    }

    plot_each_F(object$Y, object$X, object$bpwpm_params)
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

    # If bpwpm_params are passed down to the functions
    if(class(F_mat) == "bpwpm_params"){
        # If ther's already an F_mat
        if(dim(X)[1] == dim(F_mat$estimated_F)[1]){
            cat("Old F is being used")
            F_mat <- F_mat$estimated_F
        }
        # If new data is being used
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

# Methods for ploting 2D Graphs
#-------------------------------------------------------------------------------
#' Wrapper Function for 2D Input Plots
#'
#' To better understand the model, we can visualize it, however due to universal
#' limitations, plots are only availabe for X inputs of 2 dimentions. ie: only
#' two factors included on the regresion.
#'
#' @param Y A response vector of size n
#' @param X An input matrix of size n*2.
#' @param bpwpm_params an object of the class bpwpm_params of bpwpm_prediction
#' created by the functions \code{\link{posterior_params}} or
#' \code{\link{predict.bpwpm}} respectively that contains all the info about the
#' posterior parametres of the model.
#' @param n Thinness of grid for 2D projection
#' @param m Thinness of grid for 3D projection
#' @param alpha numeric - level of transparency for 2D projection
#' @param fof0 logical - if the constant function 0 is to be plotted
#'
#' @return A series of 3 plots to help ilustrate the model
#' @export
#'
plot_2D <- function(Y, X, bpwpm_params, n = 10, alpha = 0.6,
                    m = 5, fof0 = TRUE){

    # Sanitizing Inputs
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

    if(length(Y) != dim(X)[1]){
        error("Y and X have a diferent number of observations")
        geterrmessage()
    }

    # To simplify stuff
    if(class(bpwpm_params) == 'bpwpm_prediction'){
        bpwpm_params <- bpwpm_params$bpwpm_params
    }else if(class(bpwpm_params) != 'bpwpm_params'){
        error("bpwpm_params or bpwpm_prediction objects requiered to print the plots")
    }

    # Normal Data
    p <- plot_2D_data(Y,X)
    print(p)
    readline(prompt = "Press [enter] to view next plot")
    p <- plot_2D_proj(Y, X, bpwpm_params, n, alpha)
    print(p)
    readline(prompt = "Press [enter] to view next plot")
    p <- plot_3D_proj(X = X, bpwpm_params = bpwpm_params, m =  m, fof0 = fof0)
    print(p)
    readline(prompt = "Press [enter] to view next plot")
    p <- plot_2D_F(Y = Y, bpwpm_params = bpwpm_params)
    print(p)
}

#'------------------------------------------------------------------------------
#' Scatter Plot of 2D data
#'
#' Scatter plot to visualize the data and it's corresponding groups.
#'
#' @inheritParams plot_2D
#' @param f_transform logical - To clarify whether or not, the names of the axes
#' are the f transforms or the initial X's
#'
#' @return A ggplot2 scatter plot
#' @export
#'
#' @examples (Y = rbinom(100, 1, 4), X = cbind(rnorm(100), rnorm(100)))
plot_2D_data <- function(Y,X, f_transform = FALSE){

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

    if(!f_transform){
        ggplot2::ggplot(data = X, aes(x = X[, 1], y = X[ ,2], col = Y)) +
            geom_point() + xlab("X_1") + ylab("X_2")
    }else{
        ggplot2::ggplot(data = X, aes(x = X[, 1], y = X[ ,2], col = Y)) +
            geom_point() + xlab("f(X_1)") + ylab("f(X_2)")
    }
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
#' Plots the 3D representation of the projection function
#'
#' Given the set of parmeters and the input data in 2D, this function calculates
#' and plots the wireframe on a 3D linear space defined by the input matrix X.
#'
#' @inheritParams plot_2D
#'
#' @return a 3D WireFrame Lattice Plot
#' @export
#'
plot_3D_proj <- function(X, bpwpm_params, m, fof0 = TRUE){

    if(dim(X)[2] != 2){
        error("Only a 2D plot can be made. X matrix has diferent dimensions")
        geterrmessage()
    }

    mins <- apply(X, 2, min)
    maxs <- apply(X, 2, max)

    linspace <- expand.grid(X1 = seq(mins[1], maxs[1], by = 1/m),
                            X2 = seq(mins[2], maxs[2], by = 1/m))

    linspace$eta <- model_eta(new_X = linspace,
                            bpwpm_params = bpwpm_params)

    if(fof0){
        n <- dim(linspace)[1]
        linspace$g <- rep(1, times = n)

        # Building the 0 gridspace
        linspace <- rbind(linspace, data.frame(X1 = linspace[,1],
                                               X2 = linspace[,2],
                                               eta = rep(0, times = n),
                                               g = rep(0, times = n)))
    }

    lattice::wireframe(eta ~ X1 * X2,
                       data = linspace,
                       group = g,
                       drape = TRUE,
                       aspect = c(1,1),
                       main = paste("3D plot for: M = ", bpwpm_params$M,
                                    ", J = ", bpwpm_params$J,
                                    ", K = ", bpwpm_params$K,
                                    sep = ""),
                       frame.plot = FALSE,
                       colorkey = FALSE,
                       scales = list(arrows = FALSE))
    # col.groups = rgb(c(255,0,0), c(0,255,0), alpha = 70,maxColorValue = 255),
    # col.regions = colorRampPalette(c("blue", "red"))(50))
    # at = 0, col.regions = c("red", "blue"))
}

#-------------------------------------------------------------------------------
#' Plots the transformed input space F
#'
#' @inheritParams plot_2D
#'
#' @return A plot of F1(X1) and F2(X2)
#' @export
#'
#' @examples
plot_2D_F <- function(Y, bpwpm_params){
    plot_2D_data(Y = Y, X = bpwpm_params$estimated_F[,2:3],
                 f_transform = TRUE)
}



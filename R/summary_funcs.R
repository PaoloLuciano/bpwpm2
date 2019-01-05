# Generic S3 summary functions for the package bpwpm
#-------------------------------------------------------------------------------
#' Summarizes an bpwpm run
#'
#' Once the model has been run, soem relevant statistics are given for each one
#' of the parameters to evaluate its convergence and relevance.
#'
#' @param object of the class bpwpm
#' @param digits number of significant digits for the figures
#' @param ... additional parameters
#'
#' @return Prints out summarized information
#' @export
#'
summary.bpwpm <- function(object, thin = 0, burn_in = 0,
                          digits = max(3L, getOption("digits") - 3L), ...) {

    if(!('bpwpm' %in% class(object))){
        error("Object not of the class bpwpm")
        geterrmessage()
    }

    cat(object$info)

    print_info <- function(mat, title, digits){
        cat("\n",title,"\n")
        qq <- sapply(mat, quantile)
        temp <- do.call(data.frame,
                        list(min = signif(qq[1, ], digits),
                             FirstQ. = signif(qq[2, ], digits),
                             mean = signif(sapply(mat, mean), digits),
                             Median = signif(qq[3, ], digits),
                             ThirdQ. = signif(qq[4, ], digits),
                             Max = signif(qq[5, ], digits),
                             Sd. = signif(sapply(mat, sd), digits)))
        print(t(temp))
    }

    thin_object <- thin_bpwpm(object, burn_in, thin)

    print_info(mat = thin_object$beta, title = "Betas", digits = digits)
}

#-------------------------------------------------------------------------------

#' Summarizes the parameters of an bpwpm
#'
#' @param object An object of the class bpwpm_params
#' @param digits Number of significant digits for the figures
#' @param verb To print out additional info
#' @param ... Arguments to be passed to or from other methods
#'
#' @return Prints out info
#' @export
#'
summary.bpwpm_params <- function(object,
                                 digits = max(3L, getOption("digits") - 3L),
                                 verb = FALSE, ...){

    if(!('bpwpm_params' %in% class(object))){
        error("Object not of the class bpwpm_params")
        geterrmessage()
    }

    cat("\nPosterior Estimated Params: \n")

    print(t(signif(object$beta,digits = digits)))

    if(verb){s
        cat("\nNodes\n")
        print(object$tau)

        # cat("\nF_Transformation\n")
        # print(object$params$estimated_F)
    }
}

#-------------------------------------------------------------------------------

#' Summarizes an bpwpm prediciton
#'
#' @param object An object of the class bpwpm_prediction
#' @param digits Number of significant digits for the figures
#' @param verb To print out additional info
#' @param ... Arguments to be passed to or from other methods
#'
#' @return Prints out info
#' @export
#'
summary.bpwpm_prediction <- function(object,
                                     digits = max(3L, getOption("digits") - 3L),
                                     verb = FALSE, ...){

    if(!('bpwpm_prediction' %in% class(object))){
        error("Object not of the class bpwpm_prediction")
        geterrmessage()
    }

    cat(object$info)
    cat("Prediction Results:\n",
        "\tAccuracy:\t", signif(object$accuracy*100, (digits - 1)), "%\n",
        "\tLog-Loss:\t", signif(object$log_loss, digits), "\n",
        "\tType of Posterior: ", object$type, "\n",
        "\tConfusion Matrix:", "\n")
    print(object$contingency_table)

    cat("\nPosterior Estimated Params: \n")
    print(t(signif(object$bpwpm_params$beta,digits = digits)))

    if(verb){
        cat("\nNodes\n")
        print(t(signif(object$bpwpm_params$tau,digits = digits)))

        # cat("\nF_Transformation\n")
        # print(object$bpwpm_params$estimated_F)
    }
}

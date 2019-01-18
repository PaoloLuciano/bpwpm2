# Generic S3 prediction method for the bpwpm
#-------------------------------------------------------------------------------
#' Predict Method for a bpwpm
#'
#' @param object An object of the class bpwpm
#' @param new_Y New Y data to make the prediction for. However, the train Y data
#' can be used to evaluate the model
#' @param new_X New X data to make the prediction for. However, the train X data
#' can be used to evaluate the model
#' @inheritParams posterior_params
#' @param ... further arguments passed to or from other methods
#'
#' @return An object of the class \code{bpwpm_prediction} containing:
#' \describe{
#'   \item{Info}{A formated string that describes the basics of the model}
#'   \item{type}{The type of posterior probability used}
#'   \item{fitted_probabilities}{The fitted probabilities for the input}
#'   \item{bpwpm_params}{An object of the class \code{bpwpm_params}, created by
#'   the function \code{\link{posterior_params}}}
#'   \item{contingency_table}{The confusion matrix for this model, and
#'   prediction}
#'   \item{accuracy}{The accuracy of the model. This can be misleading.}
#'   \item{log_loss}{The log loss for the Y response and the fitted
#'   probabilities}
#'   \item{X}{The X input matrix passed down to the method. Used for plotting
#'   methods}
#'   \item{Y}{The Y input matrix passed down to the method. Used for plotting
#'   methods}
#' }
#' @export
#'
#' @examples (model1, train_Y, train_X, 2, 1000, mean)
#' (model2, test_Y, test_X, 1, 0, mode)

predict.bpwpm <- function(object, new_Y, new_X,
                          burn_in = 0, thin = 0, type = 'mean', ...) {

    if(!('bpwpm' %in% class(object))){
        error("Object not of the class bpwpm")
        geterrmessage()
    }

    if(class(new_Y) == "factor"){
        new_Y <- as.integer(new_Y) - 1
    }

    post_params <- posterior_params(object, burn_in = burn_in, thin = thin, type)
    p <- posterior_probs(new_X, post_params)

    model_predict <- list(info = object$info,
                          type = type,
                          bpwpm_params = post_params,
                          contingency_table = contingency_table(new_Y, p),
                          accuracy = accuracy(new_Y, p),
                          log_loss = log_loss(new_Y, p, verb = FALSE),
                          fitted_probabilities = p,
                          X = new_X,
                          Y = new_Y)

    class(model_predict) <- "bpwpm_prediction"

    return(model_predict)
}

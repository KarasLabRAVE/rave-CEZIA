#' Calculate adjacency matrices and fragility matrix from iEEG recording
#'
#' The function calculates the neural fragility column
#' from an adjacency matrix in each time window
#'
#' @source Recreation of the method described in
#' Li A, Huynh C, Fitzgerald Z, Cajigas I, Brusko D, Jagid J, et al.
#' Neural fragility as an EEG marker of the seizure onset zone.
#'  Nat Neurosci. 2021 Oct;24(10):1465–74
#' (\href{https://pubmed.ncbi.nlm.nih.gov/34354282/}{pubmed}).
#' We have found solutions to fill up missing details in the paper method description
#'
#' @param ieegts Numeric. A matrix of iEEG time series x(t),
#' with time points as rows and electrodes names as columns
#' @param window Integer. The number of time points to use in each window
#' @param step Integer. The number of time points to move the window each time
#' @param lambda Numeric. The lambda value to use in the ridge regression.
#' If NULL, the lambda will be chosen automatically
#' ensuring that ensuring that the adjacent matrix is stable (see details)
#' @param nSearch Integer. Number of minimization to compute the fragility row
#'
#' @return A list containing the normalized ieegts,
#' adjacency matrices, fragility, and R^2 values
#'
#' @examples
#' ## A simple example
#' data <- matrix(rnorm(100), nrow = 10)
#' window <- 10
#' step <- 5
#' lambda <- 0.1
#' calcAdjFrag(ieegts = data, window = window,
#' step = step, lambda = lambda)
#'
#' ## A more realistic example, but it will take a while to run
#' \dontrun{
#' data("pt01Epochm1sp2s")
#' window <- 250
#' step <- 125
#' lambda <- NULL
#' nSearch <- 100
#' title <- "PT01 seizure 1"
#' resfrag <- calcAdjFrag(ieegts = pt01Epochm1sp2s, window = window,
#'   step = step, lambda = lambda,nSearch=nSearch)
#' }
#'
#'
#' @details
#' 1/ For each time window i, a discrete stable Linear time system
#' (adjacency matrix) is computed named \eqn{A_i}
#' such that
#' \eqn{A_i x(t) = x(t+1)}
#' option Lambda=NULL ensures that the matrix is stable
#'
#' 2/For each stable estimated \eqn{A_i}, the minimum norm perturbation \eqn{\Gamma_{ik}} (k index of the electrodes)
#' for column perturbation is computed.
#' Each column is normalized \eqn{\frac{max(\Gamma_{i})-\Gamma_{ik}}{max(\Gamma_i)}}
#'
#' @export
calcAdjFrag <- function(ieegts, window, step, lambda = NULL, nSearch=100) {
    ## check the input types
    stopifnot(isWholeNumber(window))
    stopifnot(isWholeNumber(step))
    stopifnot(is.null(lambda) | is.numeric(lambda))

    ## The input matrix must have at least window rows
    stopifnot(nrow(ieegts) >= window)


    ## Number of electrodes and time points
    n_tps <- nrow(ieegts)
    n_elec <- ncol(ieegts)

    electrodeList <- colnames(ieegts)

    # Number of steps
    nSteps <- floor((n_tps - window) / step) + 1

    scaling <- 10^floor(log10(max(ieegts)))
    ieegts <- ieegts / scaling

    ## create adjacency array (array of adj matrices for each time window)
    ## iw: The index of the window we are going to calculate fragility
    res <- lapply(seq_len(nSteps), function(iw) {
        ## Sample indices for the selected window
        si <- seq_len(window - 1) + (iw - 1) * step
        ## measurements at time point t
        xt <- ieegts[si, ]
        ## measurements at time point t plus 1
        xtp1 <- ieegts[si + 1, ]

        ## Coefficient matrix A (adjacency matrix)
        ## each column is coefficients from a linear regression
        ## formula: xtp1 = xt*A + E
        if (is.null(lambda)) {
            Ai <- ridgesearchlambdadichomotomy(xt, xtp1, intercept = FALSE)
        } else {
            Ai <- ridge(xt, xtp1, intercept = FALSE, lambda = lambda)
        }

        R2 <- ridgeR2(xt, xtp1, Ai)

        list(Ai = Ai, R2 = R2)
    })

    A <- unlist(lapply(res, function(w) {
        w$Ai
    }))
    ## TODO: Why do you want to do this? very error prone
    dim(A) <- c(n_elec, n_elec, nSteps)
    dimnames(A) <- list(
        Electrode1 = electrodeList,
        Electrode2 = electrodeList,
        Step = seq_len(nSteps)
    )

    R2 <- unlist(lapply(res, function(w) {
        w$R2
    }))
    dim(R2) <- c(n_elec, nSteps)
    dimnames(R2) <- list(
        Electrode = electrodeList,
        Step = seq_len(nSteps)
    )

    if (is.null(lambda)){
        lambdas <- sapply(res, function(w) {
            attr(w$Ai, "lambdaopt")
        })
    } else {
        lambdas <- rep(lambda, length(res))
    }



    # calculate fragility
    f <- sapply(seq_len(nSteps), function(iw) {
        fragilityRow(A[, , iw],nSearch=nSearch) # Normalized minimum norm perturbation for Gammai (time window iw)
    })
    dimnames(f) <- list(
        Electrode = electrodeList,
        Step = seq_len(nSteps)
    )

    ## TODO: Is this consistent with the method in the paper?
    # ranked fragility map
    f_rank <- matrix(rank(f), nrow(f), ncol(f))
    attributes(f_rank) <- attributes(f)
    f_rank <- f_rank / max(f_rank)

    Fragility(
        ieegts = ieegts,
        adj = A,
        frag = f,
        frag_ranked = f_rank,
        R2 = R2,
        lambdas = lambdas
    )

}

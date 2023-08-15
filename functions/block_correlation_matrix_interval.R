# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------


# function to generate diagonal blocks of a correlation matrix
diagonal_block <- function(scale_size, interval_W) {
  block <- diag(1, scale_size)
  lower <- lower.tri(block)
  upper <- upper.tri(block)
  block[lower] <- runif(scale_size * (scale_size-1) / 2,
                        min = interval_W[1], max = interval_W[2])
  block[upper] <- t(block)[upper]
  block
}

# function to generate offdiagonal blocks of a correlation matrix
offdiagonal_block <- function(scale_size, interval_B) {
  matrix(runif(scale_size^2, min = interval_B[1], max = interval_B[2]),
         nrow = scale_size, ncol = scale_size)
}


#' Generate a block correlation matrix
#'
#' @param num_scales Number of scales
#' @param scale_size Number of items in a scale (equal across scales)
#' @param interval_W Interval from which unique within-scale correlations are
#' randomly drawn
#' @param interval_B Interval from which unique between-scale correlations are
#' randomly drawn
#'
#' @return A positive-definite correlation matrix
cor_mat_block <- function(num_scales, scale_size, interval_W, interval_B) {

  ## draw correlation matrix until it can be corrected to be positive-definite
  continue_while <- TRUE
  while (continue_while) {

    # initialize correlation matrix as a list of lists
    # (outer list corresponds to columns, inner list to rows of the given column)
    out <- replicate(num_scales,
                     replicate(num_scales, NULL, simplify = FALSE),
                     simplify = FALSE)
    # loop over indices of blocks in the correlation matrix, and generate those
    # blocks following the R convention of building a matrix by column
    # (k is the index of the row, l is the index of the column)
    seq_scales <- seq_len(num_scales)
    for (l in seq_scales) {
      for (k in seq_scales) {
        if (k == l) {
          out[[l]][[k]] <- diagonal_block(scale_size, interval_W)
        } else if (k > l) {
          out[[l]][[k]] <- offdiagonal_block(scale_size, interval_B)
        } else out[[l]][[k]] <- t(out[[k]][[l]])
      }
    }
    # put correlation matrix together
    out <- do.call(cbind, lapply(out, function(column) do.call(rbind, column)))

    # compute the nearest positive-definite correlation matrix
    nearPD_object <- Matrix::nearPD(out, corr = TRUE,
                                    keepDiag = TRUE,
                                    ensureSymmetry = TRUE,
                                    base.matrix = TRUE)

    # break while loop if we found a positive-definite matrix
    # (that is, algorithm converged)
    continue_while <- !nearPD_object$converged

  }

  ## extract final correlation matrix
  nearPD_object$mat

}

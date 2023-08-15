# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------


# function to generate diagonal blocks of a correlation matrix
diagonal_block <- function(scale_size, rho_W) {
  block <- matrix(rho_W, nrow = scale_size, ncol = scale_size)
  diag(block) <- 1
  block
}

# function to generate offdiagonal blocks of a correlation matrix
offdiagonal_block <- function(scale_size, rho_B) {
  matrix(rho_B, nrow = scale_size, ncol = scale_size)
}


#' Generate a block correlation matrix
#'
#' @param num_scales Number of scales
#' @param scale_size Number of items in a scale (equal across scales)
#' @param rho_W Value of the within-scale correlations
#' @param rho_B Value of the between-scale correlations
#'
#' @return A positive-definite correlation matrix
cor_mat_block <- function(num_scales, scale_size, rho_W, rho_B) {

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
          out[[l]][[k]] <- diagonal_block(scale_size, rho_W)
        } else if (k > l) {
          out[[l]][[k]] <- offdiagonal_block(scale_size, rho_B)
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

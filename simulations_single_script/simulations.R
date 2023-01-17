# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------


# load packages, but functions should still be called with package::function()
library("parallel")
library("simstudy")
library("ccaPP")


# control parameters for data generation
n <- 1000                               # number of observations
p <- 2                                  # number of items
prob <- c(0.05, 0.25, 0.4, 0.25, 0.05)  # probabilities of response categories
L <- length(prob)                       # number of response categories
rho <- 0.7                              # target correlation between items
R <- 1000                               # number of simulation runs
seed <- 20230111                        # seed of the random number generator

# define matrix with probabilities of response categories per item
prob_mat <- matrix(prob, nrow = p, ncol = L, byrow = TRUE)

# define correlation matrix
Rho <- matrix(rho, nrow = p, ncol = p)
diag(Rho) <- 1

# control parameters for random respondents
epsilons <- seq(0, 0.3, by = 0.05)
epsilon_max <- max(epsilons)


# it is very easy to use parallel computing on Unix systems, but not on Windows
if (.Platform$OS.type == "windows") {
  n_cores <- 1              # use only one CPU core
} else {
  n_cores <- 2              # number of CPU cores to be used
  RNGkind("L'Ecuyer-CMRG")  # use parallel random number streams
}


# run simulation
cat(paste(Sys.time(), ": starting ...\n"))
set.seed(seed)
results_list <- parallel::mclapply(seq_len(R), function(r) {

  # print simulation run
  cat(paste(Sys.time(), sprintf(":   run = %d\n", r)))

  # initialize data table
  initial <- simstudy::genData(n)
  # generate correlated rating-scale items
  data <- simstudy::genOrdCat(initial, baseprobs = prob_mat, prefix = "item",
                              corMatrix = Rho)
  # drop ID and convert to integer matrix
  data <- as.matrix(data[, -1])
  storage.mode(data) <- "integer"

  # generate probabilities of being a random respondents
  careless_probabilities <- runif(n)

  # order observations according to probabilities of being random respondents,
  # which makes it easier to keep previous careless respondents the same as the
  # contamination level increases (for maximum comparability)
  order <- order(careless_probabilities)
  data <- data[order, ]
  careless_probabilities <- careless_probabilities[order]

  # generate random responses to be used for careless respondents
  n_careless_max <- sum(careless_probabilities < epsilon_max)
  data_careless <- replicate(p, sample.int(L, n_careless_max, replace = TRUE))

  # loop over contamination levels
  results_r <- lapply(epsilons, function(epsilon) {

    # turn selected observations into careless respondents: since the
    # observations are sorted according to the probability of being careless,
    # this keeps previous careless respondents the same as the contamination
    # level increases
    if (epsilon > 0) {
      careless <- which(careless_probabilities < epsilon)
      data[careless, ] <- data_careless[careless, ]
    }

    # compute Pearson correlation
    df_pearson <- tryCatch({
      pearson <- ccaPP::corPearson(data[, 1], data[, 2])
      data.frame(Run = r, epsilon = epsilon, Method = "Pearson",
                 Correlation = pearson)
    }, error = function(e) NULL, warning = function(w) NULL)

    # compute Kendall correlation
    df_kendall <- tryCatch({
      kendall <- ccaPP::corKendall(data[, 1], data[, 2])
      data.frame(Run = r, epsilon = epsilon, Method = "Kendall",
                 Correlation = kendall)
    }, error = function(e) NULL, warning = function(w) NULL)

    # combine results
    rbind(df_pearson, df_kendall)

  })

  # combine results from current simulation run into data frame
  do.call(rbind, results_r)

}, mc.cores = n_cores)

# combine results into data frame
results <- do.call(rbind, results_list)

# save results to file
file_results <- "simulations_single_script/results/results_n=%d.RData"
save(results, n, p, prob, rho, seed, file = sprintf(file_results, n))

# print message that simulation is done
cat(paste(Sys.time(), ": finished.\n"))


# aggregate results over the simulation runs
library("dplyr")
aggregated <- results %>%
  group_by(epsilon, Method) %>%
  summarize(Correlation = mean(Correlation),
            .groups = "drop")

# plot average results over the simulation runs
library("ggplot2")
p <- ggplot() +
  geom_line(aes(x = epsilon, y = Correlation, color = Method),
            data = aggregated)

# save plot to file
file_plot <- "simulations_single_script/figures/results_n=%d.pdf"
pdf(file = sprintf(file_plot, n), width = 5, height = 3.5)
print(p)
dev.off()

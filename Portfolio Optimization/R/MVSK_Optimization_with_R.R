# Step 1: Import the necessary libraries

# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(corrplot)
library(moments)
library(quantmod)
library(PerformanceAnalytics)
library(DEoptim)
library(parallel)

# Define the tickers for MAANG stocks
tickers <- c("META", "AAPL", "AMZN", "NFLX", "GOOGL")

# Download the stock data from Yahoo Finance
getSymbols(tickers, src = "yahoo", from = "2015-01-01", to = "2025-01-01")

# Extract the adjusted closing prices for each stock
data <- do.call(merge, lapply(tickers, function(ticker) Cl(get(ticker))))

# Compute and display daily log returns
returns <- na.omit(diff(log(data)))
print(head(returns))


# Step 2: Manually calculate mean, variance, skewness, coskewness, kurtosis

# Mean returns
mean_returns <- colMeans(returns)

# Variance
variance <- apply(returns, 2, var)

# Skewness
skewness_vals <- apply(returns, 2, skewness)

# Kurtosis
kurtosis_vals <- apply(returns, 2, kurtosis)

# Coskewness (3rd central moment)
coskewness_matrix <- matrix(NA, nrow = ncol(returns), ncol = ncol(returns))
for (i in 1:ncol(returns)) {
  for (j in 1:ncol(returns)) {
    coskewness_matrix[i, j] <- mean((returns[, i] - mean_returns[i]) * (returns[, j] - mean_returns[j])^2)
  }
}
colnames(coskewness_matrix) <- colnames(returns)
rownames(coskewness_matrix) <- colnames(returns)

# Cokurtosis (4th central moment)
cokurtosis_matrix <- matrix(NA, nrow = ncol(returns), ncol = ncol(returns))
for (i in 1:ncol(returns)) {
  for (j in 1:ncol(returns)) {
    cokurtosis_matrix[i, j] <- mean((returns[, i] - mean_returns[i])^2 * (returns[, j] - mean_returns[j])^2)
  }
}
colnames(cokurtosis_matrix) <- colnames(returns)
rownames(cokurtosis_matrix) <- colnames(returns)


# Step 3: Reshape the matrices into long format for ggplot

# Reshape Covariance Matrix
covariance_matrix <- cov(returns)
cov_df <- melt(covariance_matrix)
colnames(cov_df) <- c("Asset1", "Asset2", "Covariance")

# Reshape Coskewness Matrix
coskewness_df <- melt(coskewness_matrix)
colnames(coskewness_df) <- c("Asset1", "Asset2", "Coskewness")

# Reshape Cokurtosis Matrix
cokurtosis_df <- melt(cokurtosis_matrix)
colnames(cokurtosis_df) <- c("Asset1", "Asset2", "Cokurtosis")


# Step 4: Plot the heatmaps using ggplot2

# 1. Covariance Matrix Heatmap
ggplot(cov_df, aes(x = Asset1, y = Asset2, fill = Covariance)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "purple", high = "blue", mid = "white") +
  theme_minimal() +
  labs(title = "Covariance Matrix", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, hjust = 1))

# 2. Coskewness Matrix Heatmap
ggplot(coskewness_df, aes(x = Asset1, y = Asset2, fill = Coskewness)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "pink", high = "red", mid = "white") +
  theme_minimal() +
  labs(title = "Coskewness Matrix", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, hjust = 1))

# 3. Cokurtosis Matrix Heatmap
ggplot(cokurtosis_df, aes(x = Asset1, y = Asset2, fill = Cokurtosis)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "blue", high = "green", mid = "white") +
  theme_minimal() +
  labs(title = "Cokurtosis Matrix", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, hjust = 1))


# Step 5: Function to calculate portfolio performance metrics

portfolio_performance <- function(weights, mean_returns, covariance_matrix, skewness_vals, kurtosis_vals) {
  # Portfolio return (mean)
  portfolio_return <- sum(weights * mean_returns)
  
  # Portfolio variance
  portfolio_variance <- t(weights) %*% covariance_matrix %*% weights
  
  # Portfolio skewness
  portfolio_skewness <- sum(weights * skewness_vals)
  
  # Portfolio kurtosis
  portfolio_kurtosis <- sum(weights * kurtosis_vals)
  
  return(c(portfolio_return, portfolio_variance, portfolio_skewness, portfolio_kurtosis))
}


# Step 6: Build the multi-objective optimization model

# Define the objective function (Polynomial Goal Programming)
objective_function <- function(weights) {
  
  # Normalize the weights so that they sum to 1
  weights <- weights / sum(weights)
  
  # Extract portfolio performance metrics
  portfolio_performance_metrics <- portfolio_performance(weights, mean_returns, covariance_matrix, skewness_vals, kurtosis_vals)
  portfolio_return <- portfolio_performance_metrics[1]
  portfolio_variance <- portfolio_performance_metrics[2]
  portfolio_skewness <- portfolio_performance_metrics[3]
  portfolio_kurtosis <- portfolio_performance_metrics[4]
  
  # Maximize mean and skewness, minimize variance and kurtosis
  objective = -portfolio_return + portfolio_variance - portfolio_skewness + portfolio_kurtosis
  
  return(objective)
}

# Constraints: Set up the DEoptim parameters: weights for each stock
lower_bounds <- rep(0, length(tickers))  # Lower bounds (no short selling)
upper_bounds <- rep(1, length(tickers))  # Upper bounds (maximum weight of 1)

# Run DEoptim
de_results <- DEoptim(objective_function, 
                      lower = lower_bounds, 
                      upper = upper_bounds, 
                      control = list(NP = 50, itermax = 500))

# Step 7: Display the results

# Get the optimized weights
optimized_weights <- de_results$optim$bestmem
optimized_weights <- optimized_weights / sum(optimized_weights)  # Normalize to sum to 1

# Display the optimized weights
names(optimized_weights) <- tickers
print(optimized_weights)

portfolio_returns <- Return.portfolio(returns, weights = optimized_weights)

# Performance summary
table.AnnualizedReturns(portfolio_returns)
charts.PerformanceSummary(portfolio_returns)






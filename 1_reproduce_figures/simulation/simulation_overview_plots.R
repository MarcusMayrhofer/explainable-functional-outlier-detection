# =============================================================================
# simulation_overview_plots.R
#
# One main-paper figure: an overview of the four outlier types. Simulates its
# own small example data, so it reads no stored results.
#
# IN : ../../functions/ (sourced)
# OUT: plots/overview.pdf
# =============================================================================

library(robustmatrix)
library(foreach)
library(doParallel)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(fda)
require(this.path)
setwd(this.path::this.dir())
source("../../functions/multivariate_functional_pca_function.R")
source("../../functions/multivariate_functional_data_functions.R")
source("../../functions/helper_functions_simulation.R")
source("../../functions/functions_for_simulations.R")




#MEAN FUNCTIONS
f1 <- function(t,lambda = 0.5) 30*t^(1 + 0.5*(1-lambda))*(1-t)^(1 + 0.5*(lambda))
f2 <- function(t, lambda = 1) 4*t + lambda*(-1)^rbinom(1, 1, 0.5)*(1.8-(0.02*pi)^(-0.5)*exp(-(t-runif(n = 1, min = 0.25, max = 0.75))^2/0.02))
f3 <- function(t, lambda = 1) 4*t + lambda*2*sin(8*(t+runif(n = 1, min = 0.25, max = 0.75))*pi)
mu1.0 = function(t, p) (replicate(n = p, expr = f1(t, lambda = 1)))
mu1.2 = function(t, p) (replicate(n = p, expr = f1(t, lambda = 0)))
mu2.0 = function(t, p) (replicate(n = p, expr = 4*t)) 
mu2.1 = function(t, p, lambda = 1) (replicate(n = p, expr = f2(t, lambda = lambda)))
mu3.1 = function(t, p) (replicate(n = p, expr = f3(t)))
mu4.1 = function(t, p) (replicate(n = p, expr = f3(t, lambda = 0.1)))

#COVARIANCE FUNCTIONS
K1.0 = function(s,t, rho = 0.3)  rho*exp((-abs(s-t))/rho) 
K2.0 = function(s,t, kappa = 5, nu = 0.5, sigma = 1) rSPDE::matern.covariance(h = abs(s-t), kappa = kappa, nu = nu, sigma = sigma)
K2.1 = function(s,t, kappa = 10, nu = 0.2, sigma = 1) rSPDE::matern.covariance(h = abs(s-t), kappa = kappa, nu = nu, sigma = sigma)

set.seed(123)
from = 0 
to = 1
n = 50
p = 3
q = 100
outlier_percentage = 0.1
quant <- qchisq(0.99, p*q)
cov_function_type = c("OU", "M")

Sigma <- cellWise::generateCorMat(d = p, corrType = "ALYZ")
Sigma_outlier <- Sigma

################################################################################################
########################### ISOLATED           #################################################
################################################################################################
K <- K1.0
K_outlier <- K1.0

eps_vec = c(0.1,0.2)
lambda_vec = c(2,5,10)/10

mu_function <- mu2.0
mu_function_outlier <- function(t, p ,lambda = lambda_vec[2]) mu2.1(t, p, lambda)

create_data <- msp(from = from, to = to, p = p, q = q, n = n, 
                   mu = mu_function, cov_row = Sigma, K = K)
mu = create_data$mu
cov_row = create_data$cov_row
cov_col = create_data$cov_col
cov_row_inv = solve(create_data$cov_row)
cov_col_inv = solve(create_data$cov_col)
n_outlier <- ceiling(n*outlier_percentage)
X <- create_data$X
X_clean <- X


create_data_outlier <- msp(from = from, to = to, p = p, q = q, n = n_outlier,
                           mu = mu_function_outlier, cov_row = Sigma_outlier, 
                           K = K_outlier, mu_rand = TRUE)
X_outlier <- create_data_outlier$X
outlier_index <- sample(1:n, n_outlier)
X[,,outlier_index] <- X_outlier


matplot((X[1,,]), type = "l", col = "gray", lty = 1)
matlines(X[1,,outlier_index], type = "l", col = "black", lwd = 2, lty = 1)

X1 <- X[1,,]
colnames(X1) <- 1:n
data_isolated <- data.frame((X1), check.names = FALSE) %>% 
  mutate("t" = create_data$eval_arg) %>%
  pivot_longer(-t) %>% 
  mutate("type" = factor(ifelse(name %in% outlier_index, "outlier", "regular")))

p1 <- ggplot(data_isolated, aes(x = t, y = value, group = name)) + 
  geom_line(color = "gray") + 
  geom_line(data = data_isolated %>% filter(type == "outlier"), color = "black", linewidth = 0.7) + 
  labs(x = TeX("$t$"), y = TeX("$X_1(t)$")) + 
  theme_classic()

################################################################################################
########################### COVARIANCE         #################################################
################################################################################################

set.seed(123)
kappa_vec <- c(7, 10, 15)
nu_vec <- c(0.1,0.2,0.5)
mu_function <- mu1.0
mu_function_outlier <- mu1.0

K <- K2.0
K_outlier <- function(s, t, kappa = kappa_vec[2], nu = nu_vec[2], sigma = 1) K2.1(s, t, kappa, nu, sigma)

create_data <- msp(from = from, to = to, p = p, q = q, n = n, 
                   mu = mu_function, cov_row = Sigma, K = K)
mu = create_data$mu
cov_row = create_data$cov_row
cov_col = create_data$cov_col
cov_row_inv = solve(create_data$cov_row)
cov_col_inv = solve(create_data$cov_col)
n_outlier <- ceiling(n*outlier_percentage)
X <- create_data$X
X_clean <- X
create_data_outlier <- msp(from = from, to = to, p = p, q = q, n = n_outlier,
                           mu = mu_function_outlier, cov_row = Sigma_outlier, 
                           K = K_outlier, mu_rand = TRUE)
X_outlier <- create_data_outlier$X
outlier_index <- sample(1:n, n_outlier)
X[,,outlier_index] <- X_outlier


X2 <- X[1,,]
colnames(X2) <- 1:n
data_covariance <- data.frame((X2), check.names = FALSE) %>% 
  mutate("t" = create_data$eval_arg) %>%
  pivot_longer(-t) %>% 
  mutate("type" = factor(ifelse(name %in% outlier_index, "outlier", "regular")))

p2 <- ggplot(data_covariance, aes(x = t, y = value, group = name)) + 
  geom_line(color = "gray") + 
  geom_line(data = data_covariance %>% filter(type == "outlier"), color = "black", linewidth = 0.7) + 
  labs(x = TeX("$t$"), y = TeX("$X_1(t)$")) + 
  theme_classic()


################################################################################################
########################### EIGEN              #################################################
################################################################################################
eigen_shifts <- c(6,9,12,15)
eigen_id <- c(1,10)

create_data <- msp(from = from, to = to, p = p, q = q, n = n, 
                   mu = mu_function, cov_row = Sigma, K = K)
mu = create_data$mu
cov_row = create_data$cov_row
cov_col = create_data$cov_col
cov_row_inv = solve(create_data$cov_row)
cov_col_inv = solve(create_data$cov_col)
n_outlier <- ceiling(n*outlier_percentage)
X <- create_data$X
X_clean <- X
X_eigen1 <- X
X_eigen10 <- X

create_data_outlier1 <- msp(from = from, to = to, p = p, q = q, n = n_outlier,
                           mu = mu_function_outlier, cov_row = Sigma_outlier, 
                           K = K_outlier, mu_rand = TRUE, eigen_id = eigen_id[1], gamma = eigen_shifts[2])
create_data_outlier10 <- msp(from = from, to = to, p = p, q = q, n = n_outlier,
                           mu = mu_function_outlier, cov_row = Sigma_outlier, 
                           K = K_outlier, mu_rand = TRUE, eigen_id = eigen_id[2], gamma = eigen_shifts[2])
X_outlier1 <- create_data_outlier1$X
X_outlier10 <- create_data_outlier10$X
outlier_index <- sample(1:n, n_outlier)
X_eigen1[,,outlier_index] <- X_outlier1
X_eigen10[,,outlier_index] <- X_outlier10

X3 <- X_eigen1[1,,]
colnames(X3) <- 1:n
data_eigen1 <- data.frame((X3), check.names = FALSE) %>% 
  mutate("t" = create_data$eval_arg) %>%
  pivot_longer(-t) %>% 
  mutate("type" = factor(ifelse(name %in% outlier_index, "outlier", "regular")))

p3 <- ggplot(data_eigen1, aes(x = t, y = value, group = name)) + 
  geom_line(color = "gray") + 
  geom_line(data = data_eigen1 %>% filter(type == "outlier"), color = "black", linewidth = 0.7) + 
  labs(x = TeX("$t$"), y = TeX("$X_1(t)$"))  + 
  theme_classic()

X4 <- X_eigen10[1,,]
colnames(X4) <- 1:n
data_eigen10 <- data.frame((X4), check.names = FALSE) %>% 
  mutate("t" = create_data$eval_arg) %>%
  pivot_longer(-t) %>% 
  mutate("type" = factor(ifelse(name %in% outlier_index, "outlier", "regular")))

p4 <- ggplot(data_eigen10, aes(x = t, y = value, group = name)) + 
  geom_line(color = "gray") + 
  geom_line(data = data_eigen10 %>% filter(type == "outlier"), color = "black", linewidth = 0.7) + 
  labs(x = TeX("$t$"), y = TeX("$X_1(t)$")) + 
  theme_classic()

data_all <- rbind(cbind(data_eigen1, outliers = "shift"),
                  cbind(data_eigen10, outliers = "shape"),
                  cbind(data_covariance, outliers = "covariance-induced"),
                  cbind(data_isolated, outliers = "isolated")) %>%
  mutate(outliers = factor(outliers, levels = c("shift", "shape", "covariance-induced", "isolated")))

p_all <- ggplot(data_all, aes(x = t, y = value, group = name)) + 
  geom_line(color = "gray") + 
  geom_line(data = data_all %>% filter(type == "outlier"), color = "black", linewidth = 0.7) + 
  labs(x = TeX("$t$"), y = TeX("$X_1(t)$")) + 
  theme_classic() + 
  facet_wrap(~outliers, nrow = 1)

ggsave(filename = "plots/overview.pdf", plot = p_all, height = 3, width = 9)

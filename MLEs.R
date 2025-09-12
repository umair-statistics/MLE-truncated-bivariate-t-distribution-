#load the required packages
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RTM)
theme_set(theme_bw(base_size = 11))
theme_update(plot.title = element_text(hjust = 0.5))

#### MAXIMUM LIKELIHOOD ESTIMATION ####
#### CASE 1: EQUAL DEGREES OF FREEDOM ####

set.seed(123)
n <- 10000
df_true <- c(9,9)
theta_true <- pi/6
cutoff <- 0.5
T.Pop <- rbvt(n = n, df = df_true, theta = theta_true)
T.trun <- T.Pop[T.Pop[,1] >= cutoff, , drop = FALSE]

# Negative log-likelihood (equal dfs)
NLL.Trunc <- function(pars, data){
  df <- pars[1]
  theta <- pars[2]
  
  # invalid parameter => large penalty & return
  if (!is.finite(df) || df < 2 || !is.finite(theta)) return(1e300)
  
  t1 <- data[,1]
  t2 <- data[,2]
  
  # safe call to dbvt
  ll <- tryCatch(
    dbvt(t1 = t1, t2 = t2, df = c(df, df), theta = theta,
         cutoff = cutoff, truncation = "right", log = TRUE),
    error = function(e) rep(-1e300, length(t1))
  )
  
  if (!all(is.finite(ll))) return(1e300)
  -sum(ll)
}

# optim with bounds
start <- c(df = 5, theta = 0.5)
mle_trun_eq <- optim(par = start,
                     fn = NLL.Trunc,
                     data = T.trun,
                     method = "L-BFGS-B",
                     lower = c(2, -pi/2 + 1e-6),
                     upper = c(200,  pi/2 - 1e-6),
                     control = list(parscale = c(5, 0.5)),
                     hessian = TRUE)
mle_trun_eq$par
# Optional SEs:
se_eq <- tryCatch(sqrt(diag(solve(mle_trun_eq$hessian))), error = function(e) rep(NA, length(mle_trun_eq$par)))
se_eq

# Grid eval (use the true/expected df for theta-profile)
theta.grid <- seq(-1.47, 1.47, by = 0.01)
t1 <- T.trun[,1]; t2 <- T.trun[,2]
ll.theta <- sapply(theta.grid, function(th) {
  sum(dbvt(t1 = t1, t2 = t2, df = df_true, theta = th,
           cutoff = cutoff, truncation = "right", log = TRUE))
})
df.theta <- data.frame(theta = theta.grid, loglik = ll.theta)

p1 <- ggplot(df.theta, aes(theta, loglik)) +
  geom_line() +
  geom_vline(xintercept = theta.grid[which.max(ll.theta)], lty = 2) +
  labs(title = expression("Case 1:" ~ n[1] == ~ n[2] == 9),
       x = expression(theta ~ "(In radians)"),
       y = "Log likelihood")

# For df profile (hold theta at true/estimated)
df.grid <- seq(1, 25, by = 0.1)
ll.df <- sapply(df.grid, function(d) {
  sum(dbvt(t1 = t1, t2 = t2, df = c(d,d), theta = theta_true,
           cutoff = cutoff, truncation = "right", log = TRUE))
})
df.df <- data.frame(df = df.grid, loglik = ll.df)

p2 <- ggplot(df.df, aes(df, loglik)) +
  geom_line() +
  geom_vline(xintercept = df.grid[which.max(ll.df)], lty = 2) +
  labs(title = expression("Case 1:" ~ n[1] == ~ n[2] == 9),
       x = "Degrees of freedom (n)",
       y = "Log likelihood")

ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c("(a)","(b)"),
                       font.label = list(size = 11, color = "black", face = "plain"))



#### CASE 2: UNEQUAL DEGREES OF FREEDOM ####

set.seed(123)
df_true2 <- c(15, 9)
theta_true2 <- pi/6
T.Pop <- rbvt(n = n, df = df_true2, theta = theta_true2)
T.trun <- T.Pop[T.Pop[,1] >= cutoff, , drop = FALSE]

# Negative log-likelihood (unequal dfs)
NLL.Trunc.ue <- function(pars, data) {
  n1 <- pars[1]; n2 <- pars[2]; theta <- pars[3]
  if (!is.finite(n1) || !is.finite(n2) || n1 < 2 || n2 < 2 || !is.finite(theta)) return(1e300)
  
  t1 <- data[,1]; t2 <- data[,2]
  ll <- tryCatch(
    dbvt(t1 = t1, t2 = t2, df = c(n1, n2), theta = theta,
         cutoff = cutoff, truncation = "right", log = TRUE),
    error = function(e) rep(-1e300, length(t1))
  )
  if (!all(is.finite(ll))) return(1e300)
  -sum(ll)
}

start_ue <- c(n1 = 20, n2 = 15, theta = 0.5)
mle_trun_ue <- optim(par = start_ue,
                     fn = NLL.Trunc.ue,
                     data = T.trun,
                     method = "L-BFGS-B",
                     lower = c(2, 2, -pi/2 + 1e-6),
                     upper = c(300, 300, pi/2 - 1e-6),
                     control = list(parscale = c(20, 15, 0.5)),
                     hessian = TRUE)
mle_trun_ue$par
se_ue <- tryCatch(sqrt(diag(solve(mle_trun_ue$hessian))), error = function(e) rep(NA, length(mle_trun_ue$par)))
se_ue

# Grid evaluations (fix the other df at the true value for each 1-d profile)
theta.grid <- seq(-1.47, 1.47, by = 0.01)
t1 <- T.trun[,1]; t2 <- T.trun[,2]

ll.theta <- sapply(theta.grid, function(th) {
  sum(dbvt(t1 = t1, t2 = t2, df = df_true2, theta = th,
           cutoff = cutoff, truncation = "right", log = TRUE))
})
df.theta <- data.frame(theta = theta.grid, loglik = ll.theta)

n1.grid <- seq(1, 30, by = 0.1)
ll.n1 <- sapply(n1.grid, function(n1) {
  # hold second df at true second df (9)
  sum(dbvt(t1 = t1, t2 = t2, df = c(n1, df_true2[2]), theta = theta_true2,
           cutoff = cutoff, truncation = "right", log = TRUE))
})
df.n1 <- data.frame(n1 = n1.grid, loglik = ll.n1)

n2.grid <- seq(1, 30, by = 0.1)
ll.n2 <- sapply(n2.grid, function(n2) {
  # hold first df at true first df (15)
  sum(dbvt(t1 = t1, t2 = t2, df = c(df_true2[1], n2), theta = theta_true2,
           cutoff = cutoff, truncation = "right", log = TRUE))
})
df.n2 <- data.frame(n2 = n2.grid, loglik = ll.n2)

g1 <- ggplot(df.theta, aes(theta, loglik)) + geom_line() +
  geom_vline(xintercept = theta.grid[which.max(ll.theta)], lty = 2) +
  labs(title = expression("Case 2:" ~ n[1] == 15 ~ "," ~ n[2] == 9),
       x = expression(theta ~ "(In radians)"), y = "Log likelihood")

g2 <- ggplot(df.n1, aes(n1, loglik)) + geom_line() +
  geom_vline(xintercept = n1.grid[which.max(ll.n1)], lty = 2) +
  labs(title = expression("Case 2:" ~ n[1] == 15 ~ "," ~ n[2] == 9),
       x = expression(n[1]), y = "Log likelihood")

g3 <- ggplot(df.n2, aes(n2, loglik)) + geom_line() +
  geom_vline(xintercept = n2.grid[which.max(ll.n2)], lty = 2) +
  labs(title = expression("Case 2:" ~ n[1] == 15 ~ "," ~ n[2] == 9),
       x = expression(n[2]), y = "Log likelihood")

ggarrange(g1, g2, g3, ncol = 3, nrow = 1, labels = c("(a)","(b)","(c)"),
                       font.label = list(size = 11, color = "black", face = "plain"))



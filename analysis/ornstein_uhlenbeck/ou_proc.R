library(TMB)
library(tidyverse)
library(here)
devtools::load_all(here('sdelib'))

set.seed(10)   # Set seed for reproducibility

# Directories -------------------------------------------------------------
fig_dir <- here('analysis', 'ornstein_uhlenbeck','fig')
dir.create(fig_dir, showWarnings = F, recursive = T)

# Parameters --------------------------------------------------------------
n <- 1000                           # Grid size
t_max <- 20                         # End Time
dt <- t_max/n                       # Time discretization
t_grid <- seq(0, t_max, by = dt)    # Time Grid
mu <- 10                            # Drift coefficient
x_t0 <- 1                           # Initial value of X(t) at time 0
theta <- 0.5                        # Rate parameter
omega <- 1                          # Wiener Process sd coefficient
sigma_y <- 0.1                      # Observation SD
i_t_x_obs <- seq(1, n, by = 5)      # Observation times


# Simulate process and observations ---------------------------------------
x <- r_ou(t_max, x_t0, dt, theta, omega, mu)
y <- x[i_t_x_obs] + rnorm(length(i_t_x_obs), sd = sigma_y)


# Simulate a realization to use as initial value
x_init <- r_ou(t_max, x_t0, dt, theta, omega, mu)


# Prepare inputs for TMB --------------------------------------------------
# Data
dat <- list(t_grid = t_grid,
            i_t_x_obs = i_t_x_obs - 1,
            y = y)

# Initial Parameter Values
param <- list(theta = theta + rnorm(1),
              log_mu = log(mu) + rnorm(1),
              log_sigma_y = log(sigma_y) + rnorm(1),
              log_omega = log(omega) + rnorm(1),
              x = x_init)

# Setup optimization Function ------------------------------------------------
obj <- TMB::MakeADFun(data = c(tmb_model = 'ou', dat),
                      parameters = param,
                      random = c("x"),
                      DLL="sdelib_TMBExports", silent =F)

# If inner problem is quadratic, set smartsearch=F
newtonOption(obj, smartsearch=FALSE)

# Optimize ----------------------------------------------------------------
# Estimate parameters and latent variables
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
# Get estimates of standard errors
sdr <- sdreport(obj)
pred <- summary(sdr,"random")

# Bundle predictions into tibble ------------------------------------------
# Confidence intervals
ci_df <- tibble(t = t_grid,
             x_pred = pred[,1],
             x_sd = pred[,2]) %>%
  mutate(lb = x_pred - qnorm(0.975)*x_sd,
         ub = x_pred + qnorm(0.975)*x_sd)
# Observations
obs_df <- tibble(t = t_grid[i_t_x_obs],
                 y = y)
# Predictions
pred_df <- tibble(t = t_grid,
             x = x,
             x_pred = pred[,1]) %>%
          pivot_longer(cols = c('x','x_pred'),
                       names_to = c("Type")) %>%
          mutate(Type = fct_recode(as.factor(Type),
                                   "Ground Truth" = "x",
                                   "Predicted" = "x_pred"))

# Plot result -------------------------------------------------------------
plt <- ggplot(pred_df) +
  geom_ribbon(data = ci_df,aes(x = t, ymin= lb, ymax = ub), alpha = 0.5) +
  geom_line(aes(x = t, y = value, color = Type), size = 0.5) +
  geom_point(data = obs_df, aes(x = t, y =y), size = 0.5) +
  geom_hline(yintercept = mu) +
  labs(x = "Time",
       y = "OU Process",
       color = "") +
  theme_bw() +
  theme(text = element_text(size = 15,face = "bold"))

ggsave(file.path(fig_dir, 'ou_example.png'), height = 9,  width = 16, scale = 1/1.5)

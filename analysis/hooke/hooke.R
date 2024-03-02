library(TMB)
library(tidyverse)
library(here)
devtools::load_all(here('sdelib'))

set.seed(10)   # Set seed for reproducibility

fig_dir <- here('analysis', 'hooke','fig')
dir.create(fig_dir, showWarnings = F, recursive = T)

# Parameters --------------------------------------------------------------
n <- 1000                           # Grid size
t_max <- 50                         # End Time
dt <- t_max/n                       # Time discretization
t_grid <- seq(0, t_max, by = dt)    # Time Grid
x_init <- c(0.1, 0.1 - dt/10)
k <- 0.1                            # y coef
b <- 0.1                            # y' coef
m <- 0.1                            # y'' coef
omega <- 0.01                       # Wiener Process sd coefficient
sigma_y <- 0.05                      # Observation SD
i_t_x_obs <- seq(1, n, by = 5)      # Observation times


# Simulate process and observations ---------------------------------------
x <- r_hooke(t_max, x_init, dt, k, b, m, omega)
y <- x[i_t_x_obs] + rnorm(length(i_t_x_obs), sd = sigma_y)
plot(t_grid, x, type = "l")
points(t_grid[i_t_x_obs], y)


# Prepare inputs for TMB --------------------------------------------------
# Data
dat <- list(t_grid = t_grid,
            i_t_x_obs = i_t_x_obs - 1,
            y = y)

# Initial Parameter Values
param <- list(x = x + rnorm(length(x), 0, sd(x)),
              log_k = log(k),
              log_b = log(b),
              log_m = log(m),
              log_sigma_y = log(sigma_y),
              log_omega = log(omega))
# Fixed parameters
map <- list(log_sigma_y = factor(NA))
            # log_k = factor(NA))# ,
            # log_m = factor(NA))
            # log_b = factor(NA))

# Setup optimization Function ------------------------------------------------
obj <- TMB::MakeADFun(data = c(tmb_model = 'hooke', dat),
                      parameters = param,
                      map = map,
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

exp(opt$par)

# Plot result -------------------------------------------------------------
plt <- ggplot(pred_df) +
  geom_ribbon(data = ci_df, aes(x = t, ymin= lb, ymax = ub), alpha = 0.5) +
  geom_line(aes(x = t, y = value, color = Type, group = Type), size = 0.5) +
  geom_point(data = obs_df, aes(x = t, y =y), size = 0.5) +
  labs(x = "Time",
       y = "Hooke SDE",
       color = "") +
  theme_bw() +
  theme(text = element_text(size = 15,face = "bold"))
plt


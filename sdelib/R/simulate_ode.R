#' Euler scheme for simulating from OU Process
#'
#' @param t_max maximum time to simulate x(t)
#' @param x_t0 x at time 0
#' @param dt delta t, time increment
#' @param theta rate parameter
#' @param omega variance of BM over 1 unit time
#' @param mu Asymptotic value (drift parameter)
#'
#' @return (vector) SDE Euler solution x(t), t = 0, dt, 2*dt, ...
#' @export
r_ou <- function(t_max, x_t0, dt, theta, omega, mu){
  n <- t_max/dt + 1
  x <- 1:n
  x[1] <- x_t0
  for(i in 1:(n-1)){
    x[i+1] <- x[i] + dt*theta*(mu-x[i]) +
      sqrt(dt)*omega*rnorm(1)
  }
  return(x)
}

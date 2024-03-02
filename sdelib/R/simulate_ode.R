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


#' Euler scheme for simulating from Hooke Process
#' @param t_max maximum time to simulate x(t)
#' @param x_init x at time 0 and time dt
#' @param dt delta t, time increment
#' @param k x coef
#' @param b x' coef
#' @param m x'' coef
#' @param omega variance of BM over t = 1
#'
#' @export
r_hooke <- function(t_max, x_init, dt, k, b, m, omega){
  if(length(x_init)!=2){
    stop("Exception: length(x_init) != 2")
  }
  n <- t_max/dt + 1
  x <- 1:n
  x[1:2] <- x_init
  for(i in 2:(n-1)){
    x[i+1] <- -(b/m)*(x[i] - x[i-1])*dt -(k/m)*x[i]*dt^2 +
      2*x[i] - x[i-1] +
     sqrt(dt)*omega*rnorm(1)
  }
  return(x)

}

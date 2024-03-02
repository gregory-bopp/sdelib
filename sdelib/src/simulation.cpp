#include "RcppArmadillo.h"

// [[Rcpp::export]]
arma::dvec r_ou(const double &t_max,
                const double &x_t0,
                const double &dt,
                const double &theta,
                const double &omega,
                const double &mu
                            ) {
  int n = std::round(t_max/dt) + 1;
  arma::dvec x(n);
  x(0) = x_t0;
  for(int i = 0; i<(n-1); i++){
    x(i+1) = x(i) + dt*theta*(mu - x(i)) +
                  (sqrt(dt)*omega*norm_rand());
  }

  return x;
}

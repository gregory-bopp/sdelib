// Inference on OU stochastic differential equation.
//
// dX = - theta(mu - X)*dt + omega*dB
//
// based on discrete observations
//
// Y(i) = X(t(i)) + e(i)
//
// where e(i) is N(0, sigma_y^2)
//
// Latent variables are the states.
//
// We use Euler approximation to evalaute transition densities. The time mesh for this
// discretization is finer than the sample interval, i.e. some (many) states are unobserved.
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
// (here it's ou)
template <class Type>
Type ou(objective_function<Type>* obj) {
  DATA_VECTOR(t_grid);        // Time points of Euler approx
  DATA_VECTOR(i_t_x_obs);     // Indeces into t_grid where X is observed
  DATA_VECTOR(y);             // Observations taken. Must have same length as i_t_x_obs.

  PARAMETER_VECTOR(x);         // States at t_grid. Length = length(t_grid)
  PARAMETER(theta);            // Rate parameter in the SDE
  PARAMETER(log_mu);
  PARAMETER(log_sigma_y);
  PARAMETER(log_omega);

  // Derived Quantities
  Type sigma_y=exp(log_sigma_y);
  Type omega=exp(log_omega);
  Type mu=exp(log_mu);
  vector<Type> sigma_x(t_grid.size());

  Type ans=0;  // ans will be the resulting likelihood

  vector<Type> dt(t_grid.size()-1);
  for(int i=0;i<dt.size();i++)
    dt(i) = t_grid(i+1)-t_grid(i);

  // Include likelihood contributions from state transitions
  for(int i=0;i<dt.size();i++){
    ans -= dnorm(x(i+1),
                 x(i) + theta*(mu - x(i))*dt(i),
                 omega*sqrt(dt(i)),
                 1);
  }

  // Include likelihood contributions from measurements
  for(int i=0;i<y.size(); ++i){
    int j=CppAD::Integer(i_t_x_obs(i));
    ans-=dnorm(y(i),x(j),sigma_y,1);
  }
  return ans;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this


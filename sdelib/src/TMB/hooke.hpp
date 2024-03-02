// Inference on Hooke's Law SDE
// ODE: mx'' + bx' + kx = 0
//
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
// (here it's ou)
template <class Type>
Type hooke(objective_function<Type>* obj) {
  DATA_VECTOR(t_grid);        // Time points of Euler approx
  DATA_VECTOR(i_t_x_obs);     // Indeces into t_grid where X is observed
  DATA_VECTOR(y);             // Observations taken. Must have same length as i_t_x_obs.

  PARAMETER_VECTOR(x);         // States at t_grid. Length = length(t_grid)
  PARAMETER(log_k);
  PARAMETER(log_b);
  PARAMETER(log_m);
  PARAMETER(log_sigma_y);
  PARAMETER(log_omega);

  // Derived Quantities
  Type k = exp(log_k);
  Type b = exp(log_b);
  Type m = exp(log_m);
  Type sigma_y=exp(log_sigma_y);
  Type omega=exp(log_omega);

  Type ans=0;  // ans will be the resulting likelihood

  vector<Type> dt(t_grid.size()-1);
  for(int i=0;i<dt.size();i++)
    dt(i) = t_grid(i+1)-t_grid(i);

  // Include likelihood contributions from state transitions
  for(int i=2;i<dt.size();i++){
    ans -= dnorm(x(i+1),
                (b/m)*(x(i) + x(i-1))*dt(i)  +(k/m)*x(i)*pow(dt(i),Type(2.0)) -
                  Type(2.0)*x(i) + x(i-1),
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


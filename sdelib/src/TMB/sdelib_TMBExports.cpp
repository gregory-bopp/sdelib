// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_sdelib_TMBExports
#include <TMB.hpp>
#include "ou.hpp"
#include "hooke.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(tmb_model);
  if(tmb_model == "NULL"){
    return 0;
  }
  if(tmb_model == "ou") {
    return ou(this);
  }
  else if(tmb_model == "hooke"){
    return(hooke(this));
  }
  else {
    error("Unknown model.");
  }
  return 0;
}

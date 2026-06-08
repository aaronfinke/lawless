/*! \file simplex-lib.h simplex optimiser library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */


#ifndef SIMPLEX_LIB
#define SIMPLEX_LIB

#include <vector>


/*! Abstract base class for zero-th order function. */
class Target_fn_order_zero {
 public:
  Target_fn_order_zero() {}
  virtual ~Target_fn_order_zero() {}

  // Return number of parameters
  virtual int num_params() const = 0;

  // Return function value given the parameter vector args
  virtual double operator() ( const std::vector<double>& args ) const = 0;
};


/*! Simplex optimiser. */
class Optimiser_simplex {
 public:
  enum TYPE { NORMAL, GRADIENT };
  // Construct with cycle controls
  Optimiser_simplex( double tolerance = 0.001, int max_cycles = 50, TYPE type =
NORMAL );

  // Run the minimisation
  // On entry:
  //  target_fn               target function, class derived from Target_fn_order_zero
  //                          this function also defiens the number of parameters
  //  args[nparam+1][nparam]  list of nparam+1 start points, each a vector of nparams
  //
  // Returns:
  //  vector of refined parameters, length nparam
  //
  std::vector<double> operator() ( const Target_fn_order_zero& target_fn, const
std::vector<std::vector<double> >& args ) const;
  void debug(const int& mode) { debug_mode = mode; }
  // return number of cycles done
  int ncycles() const {return n_cycles_;}
  // return best residual
  double bestResidual() const {return best_resid_;}

 private:
  double tolerance_;
  int max_cycles_;
  TYPE type_;
  int debug_mode;
  mutable double best_resid_;
  mutable int n_cycles_;
};


#endif

#include <individual.h>
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::XPtr<process_t> exposure_process_cpp(
    Rcpp::XPtr<CategoricalVariable> health,
    double beta,
    double N,
    double dt
){
  return Rcpp::XPtr<process_t>(
    new process_t([health, beta, N, dt](size_t t){
      // Get the size of the 'I' (Infectious) state
      size_t I = health->get_size_of("I");

      // Calculate the per-capita force of infection
      double foi = beta * I / N;

      // Get a bitset of susceptible individuals
      individual_index_t S = health->get_index_of("S");

      // Sample the susceptible individuals based on the force of infection
      // Adjusted by the time step size (dt)
      Rcpp::NumericVector size_S(1);
      size_S[0] = S.size();
      Rcpp::NumericVector rate_vector(S.size());
      rate_vector = Rcpp::pexp(size_S, foi * dt, true, false); // Using pexp from Rcpp
      bitset_sample_multi_internal(S, rate_vector.begin(), rate_vector.end());

      // Queue updates to change their state from 'S' to 'I'
      health->queue_update("E", S);
    }),
    true // Set the finalizer to true to automatically delete the pointer when it's no longer used
  );
}

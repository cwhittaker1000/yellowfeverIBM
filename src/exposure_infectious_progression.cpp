#include <individual.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::XPtr<targeted_listener_t> exposure_progression_listener_cpp_internal(
    Rcpp::XPtr<CategoricalVariable> health,
    Rcpp::XPtr<TargetedEvent> exposure_progression
) {
  return Rcpp::XPtr<targeted_listener_t>(
    new targeted_listener_t([health, exposure_progression](size_t t, const individual_index_t& target) {
      std::string infectious_state = "I";
      health->queue_update(infectious_state, target);
    }),
    true // Set the finalizer to true to automatically delete the pointer when it's no longer used
  );
}

// [[Rcpp::export]]
Rcpp::XPtr<process_t> exposure_infection_progression_process_cpp(
    Rcpp::XPtr<CategoricalVariable> health,
    Rcpp::XPtr<TargetedEvent> exposure_progression,
    double gamma_shape,
    double gamma_scale,
    double dt
){
  return Rcpp::XPtr<process_t>(
    new process_t([health, exposure_progression, gamma_shape, gamma_scale, dt](size_t t){

      // Get the size of the 'E' (Exposed) state
      individual_index_t E = health->get_index_of("E");

      // Get those who are in 'E' and have already been scheduled
      individual_index_t exposed_infectious_already_scheduled = exposure_progression->get_scheduled();

      // Subset to only those who haven't been scheduled
      E &= !exposed_infectious_already_scheduled;

      // Sample infection times and then schedule them
      Rcpp::NumericVector infection_times(exposed_infectious_already_scheduled.size());
      infection_times = Rcpp::rgamma(exposed_infectious_already_scheduled.size(), gamma_shape, gamma_scale);
      std::vector<double> infection_times_std(infection_times.begin(), infection_times.end());
      exposure_progression->schedule(E, infection_times_std);

    }),
    true // Set the finalizer to true to automatically delete the pointer when it's no longer used
  );
}

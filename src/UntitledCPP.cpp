/* a little comment above */
#include <Rcpp.h>
#include <individual.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::XPtr<process_t> multi_probability_multinomial_process_internal(
    Rcpp::XPtr<CategoricalVariable> variable,
    const std::string source_state,
    const std::vector<std::string> destination_states,
    const Rcpp::XPtr<DoubleVariable> rate_variable,
    const std::vector<double> destination_probabilities
){
  // array of cumulative probabilities
  std::vector<double> cdf(destination_probabilities);
  std::partial_sum(destination_probabilities.begin(),destination_probabilities.end(),cdf.begin(),std::plus<double>());

  // make pointer to lambda function and return XPtr to R
  return Rcpp::XPtr<process_t>(
    new process_t([variable,source_state,destination_states,rate_variable,cdf](size_t t){

      // sample leavers with their unique prob
      individual_index_t leaving_individuals(variable->get_index_of(source_state));
      std::vector<double> rate_vector = rate_variable->get_values(leaving_individuals);
      bitset_sample_multi_internal(leaving_individuals, rate_vector.begin(), rate_vector.end());

      // empty bitsets to put them (their destinations)
      std::vector<individual_index_t> destination_individuals;
      size_t n = destination_states.size();
      for (size_t i=0; i<n; i++) {
        destination_individuals.emplace_back(leaving_individuals.max_size());
      }

      // random variate for each leaver to see where they go
      const auto random = Rcpp::runif(leaving_individuals.size());
      auto random_index = 0;
      for (auto it = std::begin(leaving_individuals); it != std::end(leaving_individuals); ++it) {
        auto dest_it = std::upper_bound(cdf.begin(), cdf.end(), random[random_index]);
        int dest = std::distance(cdf.begin(), dest_it);
        destination_individuals[dest].insert(*it);
        ++random_index;
      }

      // queue state updates
      for (size_t i=0; i<n; i++) {
        variable->queue_update(destination_states[i], destination_individuals[i]);
      }

    }),
    true
  );
}

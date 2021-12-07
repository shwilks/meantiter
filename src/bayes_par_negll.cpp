
#include <RcppArmadillo.h>

// [[Rcpp::export]]
double gmt_bayes_par_negll(
    arma::vec pars,
    const arma::vec &max_titers,
    const arma::vec &min_titers,
    double prior_mean_mean,
    double prior_mean_sd,
    double prior_sd_shape,
    double prior_sd_scale
) {
  
  // Setup for calculating the negative log likelihood
  double total_negll = 0;
  
  // Calculate the prior likelihood of the data
  total_negll -= R::dnorm4(pars[0], prior_mean_mean, prior_mean_sd, 1); // Mean prior
  total_negll -= R::dgamma(1.0 / (pars[1]*pars[1]), prior_sd_shape, prior_sd_scale, 1); // SD prior
  
  // Calculate the likelihood of the given the model
  for(arma::uword i = 0; i < min_titers.n_elem; ++i) {
    
    if (max_titers[i] == min_titers[i]) {
      
      total_negll -= R::dnorm4(max_titers[i], pars[0], pars[1], 1);
      
    } else {
      
      total_negll -= R::logspace_sub(
        R::pnorm5(max_titers[i], pars[0], pars[1], 1, 1),
        R::pnorm5(min_titers[i], pars[0], pars[1], 1, 1)
      );
      
    }
    
  }
  
  // Return the negative log likelihood
  return(total_negll);
  
}
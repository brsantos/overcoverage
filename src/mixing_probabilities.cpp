#define _USE_MATH_DEFINES
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]


// Inverse logit function
static inline double inv_logit(double x) {
  return exp(x) / (1 + exp(x));
}


// Parallel worker class
struct MixingWorker : public Worker {
  const arma::mat& y;
  const arma::mat& s;
  const arma::mat& e;
  const arma::mat& i;
  const arma::mat& lambda;
  const arma::mat& epsilon1;
  const arma::mat& epsilon2;
  const arma::mat& epsilon3;
  const arma::vec& initial;
  const arma::vec& first;
  int L;
  const arma::vec& probs1;
  const arma::vec& probs2;
  const arma::mat& combins;
  const arma::vec& params;
  
  RMatrix<double> mixing;
  
  MixingWorker(const arma::mat& y,
               const arma::mat& s,
               const arma::mat& e,
               const arma::mat& i,
               const arma::mat& lambda,
               const arma::mat& epsilon1,
               const arma::mat& epsilon2,
               const arma::mat& epsilon3,
               const arma::vec& initial,
               const arma::vec& first,
               int L,
               const arma::vec& probs1,
               const arma::vec& probs2,
               const arma::mat& combins,
               const arma::vec& params,
               NumericMatrix mixing)
    : y(y), s(s), e(e), i(i), lambda(lambda), epsilon1(epsilon1), 
      epsilon2(epsilon2), epsilon3(epsilon3), initial(initial), 
      first(first), L(L), probs1(probs1), probs2(probs2),
      combins(combins), params(params), mixing(mixing) {}
  
  // Compute log-likelihood for a range of individuals
  void operator()(std::size_t begin, std::size_t end){
    
    int T = y.n_cols;
    
    arma::vec mixing_prop = {inv_logit(params(0)), 1 - inv_logit(params(0))};
    
    for (std::size_t j = begin; j < end; ++j) {
      
      // Transition matrix
      arma::mat Gamma(8, 8);
      
      // Observation matrices
      arma::mat Omega1(8, 64*(std::pow(2,L-1))+2);
      arma::mat Omega2(8, 64*(std::pow(2,L-1))+2);
      
      for (int r = 0; r < 8; ++r) {
        for (int c = 0; c < (64*(std::pow(2,L-1))+2); ++c) {
          Omega1(r, c) = 0;
          Omega2(r, c) = 0;
        }   
      }
      for (int c = 0; c < (64*(std::pow(2,L-2))); ++c) {
        Omega1(0, c) = probs1(c);
        Omega2(0, c) = probs2(c); 
      }   
      Omega1(1, (64*(std::pow(2,L-2)))) = 1;
      Omega1(2, (64*(std::pow(2,L-2))+1)) = 1;
      Omega2(1, (64*(std::pow(2,L-2)))) = 1;
      Omega2(2, (64*(std::pow(2,L-2))+1)) = 1;
      for (int c = 0; c < (64*(std::pow(2,L-2))); ++c) {
        Omega1(6, (64*(std::pow(2,L-2))+2)+c) = probs1(c);
        Omega2(6, (64*(std::pow(2,L-2))+2)+c) = probs2(c);
      }   
      
      // Assume all individuals alive and observed initially registering
      arma::vec alpha1 = initial;
      arma::vec alpha2 = initial;
      
      arma::vec weighted_likelihood(2);
      
      for (int t = (first(j)-1); t < (T-1); ++t) {
        
        int k = combins(j,t+1);
        
        Omega1(3, ((k*std::pow(2,L-2)) - 1)) = 1;
        Omega1(4, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon1(j,t+1);
        Omega1(4, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon1(j,t+1);
        Omega1(5, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon2(j,t+1);
        Omega1(5, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon2(j,t+1);
        Omega1(7, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon3(j,t+1);
        Omega1(7, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon3(j,t+1);
        
        Omega2(3, ((k*std::pow(2,L-2)) - 1)) = 1;
        Omega2(4, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon1(j,t+1);
        Omega2(4, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon1(j,t+1);
        Omega2(5, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon2(j,t+1);
        Omega2(5, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon2(j,t+1);
        Omega2(7, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon3(j,t+1);
        Omega2(7, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon3(j,t+1);
        
        Gamma(0, 0) = s(j,t) * (1 - e(j,t));
        Gamma(0, 1) = 1 - s(j,t);
        Gamma(0, 2) = lambda(j,t) * s(j,t) * e(j,t);
        Gamma(0, 3) = 0;
        Gamma(0, 4) = (1 - lambda(j,t)) * s(j,t) * e(j,t);
        Gamma(0, 5) = 0;
        Gamma(0, 6) = 0;
        Gamma(0, 7) = 0;
        Gamma(1, 0) = 0;
        Gamma(1, 1) = 0;
        Gamma(1, 2) = 0;
        Gamma(1, 3) = 0;
        Gamma(1, 4) = 0;
        Gamma(1, 5) = 0;
        Gamma(1, 6) = 0;
        Gamma(1, 7) = 1;
        Gamma(2, 0) = 0;
        Gamma(2, 1) = 0;
        Gamma(2, 2) = 0;
        Gamma(2, 3) = s(j,t) * (1 - i(j,t));
        Gamma(2, 4) = 0;
        Gamma(2, 5) = 1 - s(j,t);
        Gamma(2, 6) = s(j,t) * i(j,t);
        Gamma(2, 7) = 0;
        Gamma(3, 0) = 0;
        Gamma(3, 1) = 0;
        Gamma(3, 2) = 0;
        Gamma(3, 3) = s(j,t) * (1 - i(j,t));
        Gamma(3, 4) = 0;
        Gamma(3, 5) = 1 - s(j,t);
        Gamma(3, 6) = s(j,t) * i(j,t);
        Gamma(3, 7) = 0;
        Gamma(4, 0) = s(j,t) * i(j,t);
        Gamma(4, 1) = 0;
        Gamma(4, 2) = 0;
        Gamma(4, 3) = 0;
        Gamma(4, 4) = s(j,t) * (1 - i(j,t));
        Gamma(4, 5) = 1 - s(j,t);
        Gamma(4, 6) = 0;
        Gamma(4, 7) = 0;
        Gamma(5, 0) = 0;
        Gamma(5, 1) = 0;
        Gamma(5, 2) = 0;
        Gamma(5, 3) = 0;
        Gamma(5, 4) = 0;
        Gamma(5, 5) = 0;
        Gamma(5, 6) = 0;
        Gamma(5, 7) = 1;
        Gamma(6, 0) = s(j,t) * (1 - e(j,t));
        Gamma(6, 1) = 1 - s(j,t);
        Gamma(6, 2) = lambda(j,t) * s(j,t) * e(j,t);
        Gamma(6, 3) = 0;
        Gamma(6, 4) = (1 - lambda(j,t)) * s(j,t) * e(j,t);
        Gamma(6, 5) = 0;
        Gamma(6, 6) = 0;
        Gamma(6, 7) = 0;
        Gamma(7, 0) = 0;
        Gamma(7, 1) = 0;
        Gamma(7, 2) = 0;
        Gamma(7, 3) = 0;
        Gamma(7, 4) = 0;
        Gamma(7, 5) = 0;
        Gamma(7, 6) = 0;
        Gamma(7, 7) = 1;
        
        arma::vec temp_alpha1 = alpha1;
        alpha1 = (temp_alpha1.t() * Gamma).t();
        alpha1 %= Omega1.col(y(j,t+1)-1);
        
        arma::vec temp_alpha2 = alpha2;
        alpha2 = (temp_alpha2.t() * Gamma).t();
        alpha2 %= Omega2.col(y(j,t+1)-1);
        
      }   
      
      weighted_likelihood(0) = sum(alpha1) * mixing_prop(0);
      weighted_likelihood(1) = sum(alpha2) * mixing_prop(1);
      
      // determining mixing probabilities
      mixing(j,0) = weighted_likelihood(0) / sum(weighted_likelihood);
      mixing(j,1) = weighted_likelihood(1) / sum(weighted_likelihood);
      
    }   
  }  
}; 

// [[Rcpp::export]]
NumericMatrix mixing_probs(const arma::mat& y,
                           const arma::mat& covariates,
                           const arma::mat& age2,
                           const arma::mat& age3,
                           const arma::mat& tin2,
                           const arma::mat& tin3,
                           const arma::mat& X,
                           const arma::vec& initial,
                           const arma::vec& first,
                           int L,
                           const arma::mat& combins,
                           const arma::vec& params) {
  
  int T = y.n_cols;
  int N = y.n_rows;
  
  // Specify Delta matrix to use in survival probabilities
  arma::vec Delta1(11);
  arma::vec Delta2(11);
  arma::vec Delta3(11);
  arma::vec Delta4(11);
  arma::vec Epsilon1(11);
  arma::vec Epsilon2(11);
  arma::vec Epsilon3(11);
  for (int d = 0; d < 11; ++d) {
    Delta1(d) = params(79 + d);
    Delta2(d) = params(90 + d);
    Delta3(d) = params(101 + d);
    Delta4(d) = params(112 + d);
    Epsilon1(d) = params(123 + d);
    Epsilon2(d) = params(134 + d);
    Epsilon3(d) = params(145 + d);
  }   
  
  arma::mat s(N, T);
  arma::mat e(N, T);
  arma::mat i(N, T);
  arma::mat lambda(N, T);
  arma::mat epsilon1(N, T);
  arma::mat epsilon2(N, T);
  arma::mat epsilon3(N, T);
  
  for (int j = 0; j < N; ++j) {
    for (int t = (first(j)-1); t < T; ++t) {
      s(j,t) = inv_logit(dot(covariates.row(j), Delta1) +
        (age2(j,t) * params(156)) +
        (age3(j,t) * params(157)) +
        (tin2(j,t) * params(158)) +
        (tin3(j,t) * params(159)));
      e(j,t) = inv_logit(dot(covariates.row(j), Delta2) +
        (age2(j,t) * params(160)) +
        (age3(j,t) * params(161)) +
        (tin2(j,t) * params(162)) +
        (tin3(j,t) * params(163)));
      i(j,t) = inv_logit(dot(covariates.row(j), Delta3) +
        (age2(j,t) * params(164)) +
        (age3(j,t) * params(165)) +
        (tin2(j,t) * params(166)) +
        (tin3(j,t) * params(167)));
      lambda(j,t) = inv_logit(dot(covariates.row(j), Delta4) +
        (age2(j,t) * params(168)) +
        (age3(j,t) * params(169)) +
        (tin2(j,t) * params(170)) +
        (tin3(j,t) * params(171)));
      epsilon1(j,t) = inv_logit(dot(covariates.row(j), Epsilon1) +
        (age2(j,t) * params(172)) +
        (age3(j,t) * params(173)) +
        (tin2(j,t) * params(174)) +
        (tin3(j,t) * params(175)));
      epsilon2(j,t) = inv_logit(dot(covariates.row(j), Epsilon2) +
        (age2(j,t) * params(176)) +
        (age3(j,t) * params(177)) +
        (tin2(j,t) * params(178)) +
        (tin3(j,t) * params(179)));
      epsilon3(j,t) = inv_logit(dot(covariates.row(j), Epsilon3) +
        (age2(j,t) * params(180)) +
        (age3(j,t) * params(181)) +
        (tin2(j,t) * params(182)) +
        (tin3(j,t) * params(183)));
    }
  }
  
  // Specify Beta vector to use in observation probabilities
  arma::vec Beta1(77);
  arma::vec Beta2(77);
  for (int b = 0; b < 77; ++b) {
    Beta1(b) = params(1 + b);
  }    
  for (int b = 0; b < 4; ++b) {
    Beta2(b) = params(1 + b);
  }  
  Beta2(4) = params(78);
  for (int b = 5; b < 77; ++b) {
    Beta2(b) = params(1 + b);
  }  
  
  // Register combination observation probabilities
  arma::vec logodds1(64*std::pow(2,L-2));
  arma::vec probs1(64*std::pow(2,L-2));
  for (int idx = 0; idx < (64*(std::pow(2,L-2)) - 1); ++idx) {
    logodds1(idx) = exp(dot(X.row(idx), Beta1));
  }    
  logodds1(64*(std::pow(2,L-2)) - 1) = 1; 
  probs1 = logodds1 / sum(logodds1);
  
  arma::vec logodds2(64*std::pow(2,L-2));
  arma::vec probs2(64*std::pow(2,L-2));
  for (int idx = 0; idx < (64*(std::pow(2,L-2)) - 1); ++idx) {
    logodds2(idx) = exp(dot(X.row(idx), Beta2));
  }     
  logodds2(64*(std::pow(2,L-2)) - 1) = 1; 
  probs2 = logodds2 / sum(logodds2);
  
  // Create a thread-safe probabilities container
  NumericMatrix mixing(y.n_rows, 2);
  
  // Create the worker
  MixingWorker worker(y, s, e, i, lambda, epsilon1, epsilon2, 
                      epsilon3, initial, first, L, probs1, 
                      probs2, combins, params, mixing);
  
  // Parallel execution over the range [0, y.n_rows)
  parallelFor(0, y.n_rows, worker);
  
  // Return the mixing probabilities matrix
  return mixing;
}


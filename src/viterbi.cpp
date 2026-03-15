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
struct ViterbiWorker : public Worker {
  const arma::mat& y;
  const arma::mat& X;
  const arma::vec& initial;
  const arma::vec& first;
  int L;
  const arma::mat& combins;
  const arma::mat& mixing;
  const arma::mat& s;
  const arma::mat& e;
  const arma::mat& i; 
  const arma::mat& lambda;
  const arma::mat& epsilon1;
  const arma::mat& epsilon2;
  const arma::mat& epsilon3;
  const arma::vec& probs1;
  const arma::vec& probs2;
  
  RMatrix<double> paths;
  
  ViterbiWorker(const arma::mat& y,
                const arma::mat& X,
                const arma::vec& initial,
                const arma::vec& first,
                int L,
                const arma::mat& combins,
                const arma::mat& mixing,
                const arma::mat& s,
                const arma::mat& e,
                const arma::mat& i,
                const arma::mat& lambda,
                const arma::mat& epsilon1,
                const arma::mat& epsilon2,
                const arma::mat& epsilon3,
                const arma::vec& probs1,
                const arma::vec& probs2,
                NumericMatrix paths)
    : y(y), X(X), initial(initial), first(first), L(L), combins(combins), 
      mixing(mixing), s(s), e(e), i(i), lambda(lambda), epsilon1(epsilon1), 
      epsilon2(epsilon2), epsilon3(epsilon3), probs1(probs1), probs2(probs2), 
      paths(paths) {}
  
  // Compute optimal paths for a range of individuals
  void operator()(std::size_t begin, std::size_t end){
    
    int T = y.n_cols;
    int n_states = 8;
    int n_obs = 17;
    
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
      
      // Initialise viterbi matrix and backpointer matrix
      arma::mat viterbi_matrix(n_states, n_obs);
      arma::mat backpointer(n_states, n_obs);
      // Initialisation step
      for (int state = 0; state < n_states; ++state) {
        viterbi_matrix(state, first(j)-1) = initial(state);
        backpointer(state, first(j)-1) = 0;
      }
      
      // Recursion step
      for (int t = first(j); t < T; ++t) {
        
        int k = combins(j,t);
        
        Omega1(3, ((k*std::pow(2,L-2)) - 1)) = 1;
        Omega1(4, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon1(j,t);
        Omega1(4, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon1(j,t);
        Omega1(5, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon2(j,t);
        Omega1(5, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon2(j,t);
        Omega1(7, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon3(j,t);
        Omega1(7, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon3(j,t);
        
        Omega2(3, ((k*std::pow(2,L-2)) - 1)) = 1;
        Omega2(4, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon1(j,t);
        Omega2(4, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon1(j,t);
        Omega2(5, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon2(j,t);
        Omega2(5, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon2(j,t);
        Omega2(7, ((k*std::pow(2,L-2)) - std::pow(2,L-3) - 1)) = epsilon3(j,t);
        Omega2(7, ((k*std::pow(2,L-2)) - 1)) = 1 - epsilon3(j,t);
        
        arma::mat Omega = mixing(j,0) * Omega1 + mixing(j,1) * Omega2;
        
        Gamma(0, 0) = s(j,t-1) * (1 - e(j,t-1));
        Gamma(0, 1) = 1 - s(j,t-1);
        Gamma(0, 2) = lambda(j,t-1) * s(j,t-1) * e(j,t-1);
        Gamma(0, 3) = 0;
        Gamma(0, 4) = (1 - lambda(j,t-1)) * s(j,t-1) * e(j,t-1);
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
        Gamma(2, 3) = s(j,t-1) * (1 - i(j,t-1));
        Gamma(2, 4) = 0;
        Gamma(2, 5) = 1 - s(j,t-1);
        Gamma(2, 6) = s(j,t-1) * i(j,t-1);
        Gamma(2, 7) = 0;
        Gamma(3, 0) = 0;
        Gamma(3, 1) = 0;
        Gamma(3, 2) = 0;
        Gamma(3, 3) = s(j,t-1) * (1 - i(j,t-1));
        Gamma(3, 4) = 0;
        Gamma(3, 5) = 1 - s(j,t-1);
        Gamma(3, 6) = s(j,t-1) * i(j,t-1);
        Gamma(3, 7) = 0;
        Gamma(4, 0) = s(j,t-1) * i(j,t-1);
        Gamma(4, 1) = 0;
        Gamma(4, 2) = 0;
        Gamma(4, 3) = 0;
        Gamma(4, 4) = s(j,t-1) * (1 - i(j,t-1));
        Gamma(4, 5) = 1 - s(j,t-1);
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
        Gamma(6, 0) = s(j,t-1) * (1 - e(j,t-1));
        Gamma(6, 1) = 1 - s(j,t-1);
        Gamma(6, 2) = lambda(j,t-1) * s(j,t-1) * e(j,t-1);
        Gamma(6, 3) = 0;
        Gamma(6, 4) = (1 - lambda(j,t-1)) * s(j,t-1) * e(j,t-1);
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
        
        for (int state = 0; state < n_states; ++state) {
          double max_prob = R_NegInf;
          int max_state = 0;
          for (int s_prev = 0; s_prev < n_states; ++s_prev) {
            double prob = viterbi_matrix(s_prev, t-1) * Gamma(s_prev, state) * Omega(state, y(j,t)-1);
            if (prob > max_prob) {
              max_prob = prob;
              max_state = s_prev;
            }
          }
          viterbi_matrix(state,t) = max_prob;
          backpointer(state,t) = max_state;
        }
      }
      
      // Termination step
      int best_last_state = index_max(viterbi_matrix.col(n_obs-1));
      
      // Backtracking to find the best path
      arma::vec best_path(n_obs);
      best_path(n_obs-1) = best_last_state;
      for (int t = (n_obs-2); t >= (first(j)-1); --t) {
        best_path(t) = backpointer(best_path(t+1), t+1);
      } 
      
      // log-likelihood
      for (int t = (first(j)-1); t < T; ++t) {
        paths(j,t) = best_path(t) + 1;
      }
      
    }   
  }      
};   

// [[Rcpp::export]]
NumericMatrix viterbi(const arma::mat& y,
                      const arma::mat& covariates,
                      const arma::mat& age2,
                      const arma::mat& age3,
                      const arma::mat& age4,
                      const arma::mat& tin2,
                      const arma::mat& tin3,
                      const arma::mat& combins,
                      const arma::mat& X,
                      const arma::vec& initial,
                      const arma::vec& first,
                      int L,
                      const arma::vec& params, 
                      const arma::mat& mixing) {
  
  int T = y.n_cols;
  int N = y.n_rows;
  
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
  
  arma::vec Delta1(11);
  arma::vec Delta2(11);
  arma::vec Delta3(11);
  arma::vec Delta4(11);
  arma::vec Epsilon1(11);
  arma::vec Epsilon2(11);
  arma::vec Epsilon3(11);
  for (int b = 0; b < 11; ++b) {
    Delta1(b) = params(79 + b);
    Delta2(b) = params(90 + b);
    Delta3(b) = params(101 + b);
    Delta4(b) = params(112 + b);
    Epsilon1(b) = params(123 + b);
    Epsilon2(b) = params(134 + b);
    Epsilon3(b) = params(145 + b);
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
        (age4(j,t) * params(158)) +
        (tin2(j,t) * params(159)) +
        (tin3(j,t) * params(160)));
      e(j,t) = inv_logit(dot(covariates.row(j), Delta2) +
        (age2(j,t) * params(161)) +
        (age3(j,t) * params(162)) + 
        (age4(j,t) * params(163)) +
        (tin2(j,t) * params(164)) + 
        (tin3(j,t) * params(165)));
      i(j,t) = inv_logit(dot(covariates.row(j), Delta3) +
        (age2(j,t) * params(166)) +
        (age3(j,t) * params(167)) + 
        (age4(j,t) * params(168)) +
        (tin2(j,t) * params(169)) + 
        (tin3(j,t) * params(170)));
      lambda(j,t) = inv_logit(dot(covariates.row(j), Delta4) +
        (age2(j,t) * params(171)) +
        (age3(j,t) * params(172)) + 
        (age4(j,t) * params(173)) +
        (tin2(j,t) * params(174)) + 
        (tin3(j,t) * params(175)));
      epsilon1(j,t) = inv_logit(dot(covariates.row(j), Epsilon1) +
        (age2(j,t) * params(176)) +
        (age3(j,t) * params(177)) + 
        (age4(j,t) * params(178)) +
        (tin2(j,t) * params(179)) + 
        (tin3(j,t) * params(180)));
      epsilon2(j,t) = inv_logit(dot(covariates.row(j), Epsilon2) +
        (age2(j,t) * params(181)) +
        (age3(j,t) * params(182)) + 
        (age4(j,t) * params(183)) +
        (tin2(j,t) * params(184)) + 
        (tin3(j,t) * params(185)));
      epsilon3(j,t) = inv_logit(dot(covariates.row(j), Epsilon3) +
        (age2(j,t) * params(186)) +
        (age3(j,t) * params(187)) + 
        (age4(j,t) * params(188)) +
        (tin2(j,t) * params(189)) + 
        (tin3(j,t) * params(190)));
    }
  } 
  
  // Create a thread-safe paths container
  NumericMatrix paths(N,T);
  
  // Create the worker
  ViterbiWorker worker(y, X, initial, first, L, combins, mixing, s, e, i, 
                       lambda, epsilon1, epsilon2, epsilon3, probs1, probs2, 
                       paths);
  
  // Parallel execution over the range [0, N)
  parallelFor(0, N, worker);
  
  // Sum and return the total log-likelihood
  return paths;
} 



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

// Resampling
arma::vec rmulti_rcpp(unsigned int &size,
                      NumericVector &probs) {
  int b = probs.length();
  IntegerVector outcome(b);
  rmultinom(size, probs.begin(), b, outcome.begin());
  
  return as<arma::vec>(outcome);
}

// Parallel worker class
struct LogLikelihoodWorker : public Worker {
  const arma::mat& y;
  const arma::mat& X;
  const arma::vec& initial;
  const arma::vec& first;
  int L;
  const arma::mat& combins;
  const arma::vec& mixing_prop;
  const arma::mat& s;
  const arma::mat& e;
  const arma::mat& i; 
  const arma::mat& lambda;
  const arma::mat& epsilon1;
  const arma::mat& epsilon2;
  const arma::mat& epsilon3;
  const arma::vec& probs1;
  const arma::vec& probs2;
  const arma::vec& counts;
  
  RVector<double> logL;
  
  LogLikelihoodWorker(const arma::mat& y,
                      const arma::mat& X,
                      const arma::vec& initial,
                      const arma::vec& first,
                      int L,
                      const arma::mat& combins,
                      const arma::vec& mixing_prop,
                      const arma::mat& s,
                      const arma::mat& e,
                      const arma::mat& i,
                      const arma::mat& lambda,
                      const arma::mat& epsilon1,
                      const arma::mat& epsilon2,
                      const arma::mat& epsilon3,
                      const arma::vec& probs1,
                      const arma::vec& probs2,
                      const arma::vec& counts,
                      NumericVector logL)
    : y(y), X(X), initial(initial), first(first), L(L), combins(combins), 
      mixing_prop(mixing_prop), s(s), e(e), i(i), lambda(lambda), 
      epsilon1(epsilon1), epsilon2(epsilon2), epsilon3(epsilon3), 
      probs1(probs1), probs2(probs2), counts(counts), logL(logL) {}
  
  // Compute log-likelihood for a range of individuals
  void operator()(std::size_t begin, std::size_t end){
    
    int T = y.n_cols;
    int reg_block = (L >= 2) ? (1 << (L-2)) : 1;
    int half_block = (L >= 3) ? (1 << (L-3)) : 0;
    int obs_block = 64 * reg_block;
    int obs_offset_6 = obs_block + 2;
    
    for (std::size_t j = begin; j < end; ++j) {
      double alpha1[8];
      double alpha2[8];

      for (int k_state = 0; k_state < 8; ++k_state) {
        alpha1[k_state] = initial(k_state);
        alpha2[k_state] = initial(k_state);
      }

      for (int t = (first(j)-1); t < (T-1); ++t) {

        double sj = s(j,t);
        double ej = e(j,t);
        double ij = i(j,t);
        double lj = lambda(j,t);

        double next1[8] = {0};
        double next2[8] = {0};

        // state 0 transitions
        next1[0] += alpha1[0] * sj * (1 - ej);
        next1[1] += alpha1[0] * (1 - sj);
        next1[2] += alpha1[0] * lj * sj * ej;
        next1[4] += alpha1[0] * (1 - lj) * sj * ej;

        next2[0] += alpha2[0] * sj * (1 - ej);
        next2[1] += alpha2[0] * (1 - sj);
        next2[2] += alpha2[0] * lj * sj * ej;
        next2[4] += alpha2[0] * (1 - lj) * sj * ej;

        // state 1
        next1[7] += alpha1[1];
        next2[7] += alpha2[1];

        // state 2
        next1[3] += alpha1[2] * sj * (1 - ij);
        next1[5] += alpha1[2] * (1 - sj);
        next1[6] += alpha1[2] * sj * ij;

        next2[3] += alpha2[2] * sj * (1 - ij);
        next2[5] += alpha2[2] * (1 - sj);
        next2[6] += alpha2[2] * sj * ij;

        // state 3
        next1[3] += alpha1[3] * sj * (1 - ij);
        next1[5] += alpha1[3] * (1 - sj);
        next1[6] += alpha1[3] * sj * ij;

        next2[3] += alpha2[3] * sj * (1 - ij);
        next2[5] += alpha2[3] * (1 - sj);
        next2[6] += alpha2[3] * sj * ij;

        // state 4
        next1[0] += alpha1[4] * sj * ij;
        next1[4] += alpha1[4] * sj * (1 - ij);
        next1[5] += alpha1[4] * (1 - sj);

        next2[0] += alpha2[4] * sj * ij;
        next2[4] += alpha2[4] * sj * (1 - ij);
        next2[5] += alpha2[4] * (1 - sj);

        // state 5
        next1[7] += alpha1[5];
        next2[7] += alpha2[5];

        // state 6
        next1[0] += alpha1[6] * sj * (1 - ej);
        next1[1] += alpha1[6] * (1 - sj);
        next1[2] += alpha1[6] * lj * sj * ej;
        next1[4] += alpha1[6] * (1 - lj) * sj * ej;

        next2[0] += alpha2[6] * sj * (1 - ej);
        next2[1] += alpha2[6] * (1 - sj);
        next2[2] += alpha2[6] * lj * sj * ej;
        next2[4] += alpha2[6] * (1 - lj) * sj * ej;

        // state 7
        next1[7] += alpha1[7];
        next2[7] += alpha2[7];

        for (int k_state = 0; k_state < 8; ++k_state) {
          alpha1[k_state] = next1[k_state];
          alpha2[k_state] = next2[k_state];
        }

        int obs = static_cast<int>(y(j,t+1)) - 1;
        int k = static_cast<int>(combins(j,t+1));
        int k_index = k * reg_block - 1;
        int k_false = k * reg_block - half_block - 1;

        double eps1 = epsilon1(j,t+1);
        double eps2 = epsilon2(j,t+1);
        double eps3 = epsilon3(j,t+1);

        double emit1[8] = {0};
        double emit2[8] = {0};

        if (obs >= 0) {
          if (obs < obs_block) {
            emit1[0] = probs1(obs);
            emit2[0] = probs2(obs);
          } else if (obs == obs_block) {
            emit1[1] = 1.0;
            emit2[1] = 1.0;
          } else if (obs == obs_block + 1) {
            emit1[2] = 1.0;
            emit2[2] = 1.0;
          } else if (obs >= obs_offset_6 && obs < obs_offset_6 + obs_block) {
            int idx = obs - obs_offset_6;
            emit1[6] = probs1(idx);
            emit2[6] = probs2(idx);
          }

          if (obs == k_index) {
            emit1[3] = 1.0;
            emit2[3] = 1.0;
            emit1[4] = 1.0 - eps1;
            emit2[4] = 1.0 - eps1;
            emit1[5] = 1.0 - eps2;
            emit2[5] = 1.0 - eps2;
            emit1[7] = 1.0 - eps3;
            emit2[7] = 1.0 - eps3;
          }
          if (half_block > 0 && obs == k_false) {
            emit1[4] = eps1;
            emit2[4] = eps1;
            emit1[5] = eps2;
            emit2[5] = eps2;
            emit1[7] = eps3;
            emit2[7] = eps3;
          }
        }

        for (int k_state = 0; k_state < 8; ++k_state) {
          alpha1[k_state] *= emit1[k_state];
          alpha2[k_state] *= emit2[k_state];
        }
      }

      double L1 = 0.0;
      double L2 = 0.0;
      for (int k_state = 0; k_state < 8; ++k_state) {
        L1 += alpha1[k_state];
        L2 += alpha2[k_state];
      }
      
      // log-likelihood
      double Li = mixing_prop(0) * L1 + mixing_prop(1) * L2;
      logL[j] = counts(j) * log(Li);
      
    }   
  }     
};

// --------------------------------------------------
// Legacy (matrix-based) worker for benchmarking
// --------------------------------------------------

struct LogLikelihoodWorkerLegacy : public Worker {
  const arma::mat& y;
  const arma::mat& X;
  const arma::vec& initial;
  const arma::vec& first;
  int L;
  const arma::mat& combins;
  const arma::vec& mixing_prop;
  const arma::mat& s;
  const arma::mat& e;
  const arma::mat& i; 
  const arma::mat& lambda;
  const arma::mat& epsilon1;
  const arma::mat& epsilon2;
  const arma::mat& epsilon3;
  const arma::vec& probs1;
  const arma::vec& probs2;
  const arma::vec& counts;

  RVector<double> logL;

  LogLikelihoodWorkerLegacy(const arma::mat& y,
                            const arma::mat& X,
                            const arma::vec& initial,
                            const arma::vec& first,
                            int L,
                            const arma::mat& combins,
                            const arma::vec& mixing_prop,
                            const arma::mat& s,
                            const arma::mat& e,
                            const arma::mat& i,
                            const arma::mat& lambda,
                            const arma::mat& epsilon1,
                            const arma::mat& epsilon2,
                            const arma::mat& epsilon3,
                            const arma::vec& probs1,
                            const arma::vec& probs2,
                            const arma::vec& counts,
                            NumericVector logL)
    : y(y), X(X), initial(initial), first(first), L(L), combins(combins),
      mixing_prop(mixing_prop), s(s), e(e), i(i), lambda(lambda),
      epsilon1(epsilon1), epsilon2(epsilon2), epsilon3(epsilon3),
      probs1(probs1), probs2(probs2), counts(counts), logL(logL) {}

  void operator()(std::size_t begin, std::size_t end){

    int T = y.n_cols;

    for (std::size_t j = begin; j < end; ++j) {

      arma::mat Gamma(8, 8);

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

      arma::vec individual_L(2);
      arma::vec alpha1 = initial;
      arma::vec alpha2 = initial;

      int prev_k_index = -1;
      int prev_k_false = -1;

      for (int t = (first(j)-1); t < (T-1); ++t) {

        int k = combins(j,t+1);
        int k_index = static_cast<int>((k * std::pow(2, L - 2)) - 1);
        int k_false = static_cast<int>((k * std::pow(2, L - 2)) - std::pow(2, L - 3) - 1);

        if (prev_k_index >= 0 && (prev_k_index != k_index || prev_k_false != k_false)) {
          Omega1(3, prev_k_index) = 0;
          Omega1(4, prev_k_index) = 0;
          Omega1(5, prev_k_index) = 0;
          Omega1(7, prev_k_index) = 0;
          Omega2(3, prev_k_index) = 0;
          Omega2(4, prev_k_index) = 0;
          Omega2(5, prev_k_index) = 0;
          Omega2(7, prev_k_index) = 0;

          if (prev_k_false >= 0 && prev_k_false != prev_k_index) {
            Omega1(4, prev_k_false) = 0;
            Omega1(5, prev_k_false) = 0;
            Omega1(7, prev_k_false) = 0;
            Omega2(4, prev_k_false) = 0;
            Omega2(5, prev_k_false) = 0;
            Omega2(7, prev_k_false) = 0;
          }
        }

        Omega1(3, k_index) = 1;
        Omega1(4, k_false) = epsilon1(j,t+1);
        Omega1(4, k_index) = 1 - epsilon1(j,t+1);
        Omega1(5, k_false) = epsilon2(j,t+1);
        Omega1(5, k_index) = 1 - epsilon2(j,t+1);
        Omega1(7, k_false) = epsilon3(j,t+1);
        Omega1(7, k_index) = 1 - epsilon3(j,t+1);

        Omega2(3, k_index) = 1;
        Omega2(4, k_false) = epsilon1(j,t+1);
        Omega2(4, k_index) = 1 - epsilon1(j,t+1);
        Omega2(5, k_false) = epsilon2(j,t+1);
        Omega2(5, k_index) = 1 - epsilon2(j,t+1);
        Omega2(7, k_false) = epsilon3(j,t+1);
        Omega2(7, k_index) = 1 - epsilon3(j,t+1);

        prev_k_index = k_index;
        prev_k_false = k_false;

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

      individual_L(0) = sum(alpha1) * mixing_prop(0);
      individual_L(1) = sum(alpha2) * mixing_prop(1);

      logL[j] = counts(j) * log(sum(individual_L));
    }
  }
};

// [[Rcpp::export]]
double loglikelihood_BLB_parallel(const arma::mat& y,
                              const arma::mat& covariates,
                              const arma::mat& age2,
                              const arma::mat& age3,
                              const arma::mat& age4,
                              const arma::mat& tin2,
                              const arma::mat& tin3,
                              const arma::mat& X,
                              const arma::vec& initial,
                              const arma::vec& first,
                              int L,
                              const arma::mat& combins,
                              unsigned int &n,
                              NumericVector &proportions,
                              const arma::vec& params) {
  
  int T = y.n_cols;
  int N = y.n_rows;
  
  arma::vec mixing_prop = {inv_logit(params(0)), 1 - inv_logit(params(0))};
  
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
  
  // Create resample via multinomial
  arma::vec counts = rmulti_rcpp(n, proportions);
  
  // Create a thread-safe log-likelihood container
  NumericVector logL(N);
  
  // Create the worker
  LogLikelihoodWorker worker(y, X, initial, first, L, combins, mixing_prop, 
                             s, e, i, lambda, epsilon1, epsilon2, epsilon3,
                             probs1, probs2, counts, logL);
  
  // Parallel execution over the range [0, N)
  parallelFor(0, N, worker);
  
  // Sum and return the total log-likelihood
  return -sum(logL);
} 

// [[Rcpp::export]]
double loglikelihood_BLB_parallel_legacy(const arma::mat& y,
                                     const arma::mat& covariates,
                                     const arma::mat& age2,
                                     const arma::mat& age3,
                                     const arma::mat& age4,
                                     const arma::mat& tin2,
                                     const arma::mat& tin3,
                                     const arma::mat& X,
                                     const arma::vec& initial,
                                     const arma::vec& first,
                                     int L,
                                     const arma::mat& combins,
                                     unsigned int &n,
                                     NumericVector &proportions,
                                     const arma::vec& params) {

  int T = y.n_cols;
  int N = y.n_rows;

  arma::vec mixing_prop = {inv_logit(params(0)), 1 - inv_logit(params(0))};

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

  arma::vec counts = rmulti_rcpp(n, proportions);

  NumericVector logL(N);
  LogLikelihoodWorkerLegacy worker(y, X, initial, first, L, combins, mixing_prop,
                                   s, e, i, lambda, epsilon1, epsilon2, epsilon3,
                                   probs1, probs2, counts, logL);
  parallelFor(0, N, worker);

  return -sum(logL);
}



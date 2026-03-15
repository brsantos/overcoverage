// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

static const int S = 8;

//--------------------------------------------------
// Worker
//--------------------------------------------------

struct LikelihoodWorkerCpp : public Worker {

  const arma::mat& y;
  const arma::vec& first;
  const arma::mat& combins;

  const arma::mat& s;
  const arma::mat& e;
  const arma::mat& i;
  const arma::mat& lambda;

  const arma::mat& epsilon1;
  const arma::mat& epsilon2;
  const arma::mat& epsilon3;

  const arma::vec& mixing_prop;
  const arma::vec& probs1;
  const arma::vec& probs2;
  const arma::vec& counts;

  int L;
  int T;

  RVector<double> logL;

  LikelihoodWorkerCpp(
    const arma::mat& y,
    const arma::vec& first,
    const arma::mat& combins,
    const arma::mat& s,
    const arma::mat& e,
    const arma::mat& i,
    const arma::mat& lambda,
    const arma::mat& epsilon1,
    const arma::mat& epsilon2,
    const arma::mat& epsilon3,
    const arma::vec& mixing_prop,
    const arma::vec& probs1,
    const arma::vec& probs2,
    const arma::vec& counts,
    int L,
    NumericVector logL
  )
    : y(y), first(first), combins(combins),
      s(s), e(e), i(i), lambda(lambda),
      epsilon1(epsilon1), epsilon2(epsilon2), epsilon3(epsilon3),
      mixing_prop(mixing_prop), probs1(probs1), probs2(probs2), counts(counts),
      L(L), T(y.n_cols),
      logL(logL) {}

  void operator()(std::size_t begin, std::size_t end){

    int reg_block = (L >= 2) ? (1 << (L-2)) : 1;
    int half_block = (L >= 3) ? (1 << (L-3)) : 0;
    int obs_block = 64 * reg_block;
    int obs_offset_6 = obs_block + 2;

    for(size_t j=begin; j<end; j++){

      arma::vec alpha1(S);
      arma::vec alpha2(S);

      alpha1.zeros();
      alpha2.zeros();

      alpha1(0)=1;
      alpha2(0)=1;

      for(int t=first(j)-1; t<T-1; t++){

        arma::mat Gamma(S,S,fill::zeros);

        double sj = s(j,t);
        double ej = e(j,t);
        double ij = i(j,t);
        double lj = lambda(j,t);

        // transition matrix
        Gamma(0,0)=sj*(1-ej);
        Gamma(0,1)=1-sj;
        Gamma(0,2)=lj*sj*ej;
        Gamma(0,4)=(1-lj)*sj*ej;

        Gamma(1,7)=1;

        Gamma(2,3)=sj*(1-ij);
        Gamma(2,5)=1-sj;
        Gamma(2,6)=sj*ij;

        Gamma(3,3)=sj*(1-ij);
        Gamma(3,5)=1-sj;
        Gamma(3,6)=sj*ij;

        Gamma(4,0)=sj*ij;
        Gamma(4,4)=sj*(1-ij);
        Gamma(4,5)=1-sj;

        Gamma(5,7)=1;

        Gamma(6,0)=sj*(1-ej);
        Gamma(6,1)=1-sj;
        Gamma(6,2)=lj*sj*ej;
        Gamma(6,4)=(1-lj)*sj*ej;

        Gamma(7,7)=1;

        alpha1 = (alpha1.t()*Gamma).t();
        alpha2 = (alpha2.t()*Gamma).t();

        int obs = static_cast<int>(y(j,t+1)) - 1;
        int k = static_cast<int>(combins(j,t+1));
        int k_index = k * reg_block - 1;
        int k_false = k * reg_block - half_block - 1;

        double eps1 = epsilon1(j,t+1);
        double eps2 = epsilon2(j,t+1);
        double eps3 = epsilon3(j,t+1);

        double emit1[S] = {0};
        double emit2[S] = {0};

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

        for (int k_state = 0; k_state < S; k_state++) {
          alpha1(k_state) *= emit1[k_state];
          alpha2(k_state) *= emit2[k_state];
        }

      }

      double Li =
        mixing_prop(0)*sum(alpha1) +
        mixing_prop(1)*sum(alpha2);

      logL[j] = counts(j)*std::log(Li);

    }

  }

};

// [[Rcpp::export]]
double loglikelihood_cpp(
    const arma::mat& y,
    const arma::vec& first,
    const arma::mat& combins,
    const arma::mat& s,
    const arma::mat& e,
    const arma::mat& i,
    const arma::mat& lambda,
    const arma::mat& epsilon1,
    const arma::mat& epsilon2,
    const arma::mat& epsilon3,
    const arma::vec& mixing_prop,
    const arma::vec& probs1,
    const arma::vec& probs2,
    const arma::vec& counts,
    int L
){

  int N = y.n_rows;

  NumericVector logL(N);

  LikelihoodWorkerCpp worker(
      y,first,combins,
      s,e,i,lambda,
      epsilon1,epsilon2,epsilon3,
      mixing_prop,probs1,probs2,counts,
      L,logL);

  parallelFor(0,N,worker);

  return sum(logL);

}

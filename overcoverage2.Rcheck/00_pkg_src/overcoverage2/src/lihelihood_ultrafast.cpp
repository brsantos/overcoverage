// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

static const int S = 8;

inline double inv_logit(double x){
  return 1.0/(1.0+std::exp(-x));
}


//--------------------------------------------------
// Worker
//--------------------------------------------------

struct UltraFastWorker : public Worker {

  const RMatrix<double> y;
  const RVector<double> first;
  const RMatrix<double> combins;

  const RMatrix<double> s;
  const RMatrix<double> e;
  const RMatrix<double> i;
  const RMatrix<double> lambda;

  const RMatrix<double> epsilon1;
  const RMatrix<double> epsilon2;
  const RMatrix<double> epsilon3;

  const RVector<double> mixing_prop;
  const RVector<double> probs1;
  const RVector<double> probs2;
  const RVector<double> counts;

  int L;
  int T;

  RVector<double> logL;

  UltraFastWorker(
    const NumericMatrix& y,
    const NumericVector& first,
    const NumericMatrix& combins,
    const NumericMatrix& s,
    const NumericMatrix& e,
    const NumericMatrix& i,
    const NumericMatrix& lambda,
    const NumericMatrix& epsilon1,
    const NumericMatrix& epsilon2,
    const NumericMatrix& epsilon3,
    const NumericVector& mixing_prop,
    const NumericVector& probs1,
    const NumericVector& probs2,
    const NumericVector& counts,
    int L,
    NumericVector logL
  )
    : y(y), first(first), combins(combins),
      s(s), e(e), i(i), lambda(lambda),
      epsilon1(epsilon1), epsilon2(epsilon2), epsilon3(epsilon3),
      mixing_prop(mixing_prop), probs1(probs1), probs2(probs2),
      counts(counts),
      L(L),
      T(y.ncol()),
      logL(logL) {}



  void operator()(std::size_t begin, std::size_t end){

    int reg_block = (L >= 2) ? (1 << (L-2)) : 1;
    int half_block = (L >= 3) ? (1 << (L-3)) : 0;
    int obs_block = 64 * reg_block;
    int obs_offset_6 = obs_block + 2;

    for(size_t j = begin; j < end; j++){

      double alpha1[S];
      double alpha2[S];

      for(int k=0;k<S;k++){
        alpha1[k]=0.0;
        alpha2[k]=0.0;
      }

      alpha1[0]=1.0;
      alpha2[0]=1.0;

      for(int t=first[j]-1; t<T-1; t++){

        double sj = s(j,t);
        double ej = e(j,t);
        double ij = i(j,t);
        double lj = lambda(j,t);

        double next1[S]={0};
        double next2[S]={0};


        // ---- state 0 transitions
        next1[0]+=alpha1[0]*sj*(1-ej);
        next1[1]+=alpha1[0]*(1-sj);
        next1[2]+=alpha1[0]*lj*sj*ej;
        next1[4]+=alpha1[0]*(1-lj)*sj*ej;

        next2[0]+=alpha2[0]*sj*(1-ej);
        next2[1]+=alpha2[0]*(1-sj);
        next2[2]+=alpha2[0]*lj*sj*ej;
        next2[4]+=alpha2[0]*(1-lj)*sj*ej;


        // state 1
        next1[7]+=alpha1[1];
        next2[7]+=alpha2[1];


        // state 2
        next1[3]+=alpha1[2]*sj*(1-ij);
        next1[5]+=alpha1[2]*(1-sj);
        next1[6]+=alpha1[2]*sj*ij;

        next2[3]+=alpha2[2]*sj*(1-ij);
        next2[5]+=alpha2[2]*(1-sj);
        next2[6]+=alpha2[2]*sj*ij;


        // state 3
        next1[3]+=alpha1[3]*sj*(1-ij);
        next1[5]+=alpha1[3]*(1-sj);
        next1[6]+=alpha1[3]*sj*ij;

        next2[3]+=alpha2[3]*sj*(1-ij);
        next2[5]+=alpha2[3]*(1-sj);
        next2[6]+=alpha2[3]*sj*ij;


        // state 4
        next1[0]+=alpha1[4]*sj*ij;
        next1[4]+=alpha1[4]*sj*(1-ij);
        next1[5]+=alpha1[4]*(1-sj);

        next2[0]+=alpha2[4]*sj*ij;
        next2[4]+=alpha2[4]*sj*(1-ij);
        next2[5]+=alpha2[4]*(1-sj);


        // state 5
        next1[7]+=alpha1[5];
        next2[7]+=alpha2[5];


        // state 6
        next1[0]+=alpha1[6]*sj*(1-ej);
        next1[1]+=alpha1[6]*(1-sj);
        next1[2]+=alpha1[6]*lj*sj*ej;
        next1[4]+=alpha1[6]*(1-lj)*sj*ej;

        next2[0]+=alpha2[6]*sj*(1-ej);
        next2[1]+=alpha2[6]*(1-sj);
        next2[2]+=alpha2[6]*lj*sj*ej;
        next2[4]+=alpha2[6]*(1-lj)*sj*ej;


        // state 7
        next1[7]+=alpha1[7];
        next2[7]+=alpha2[7];


        for(int k=0;k<S;k++){
          alpha1[k]=next1[k];
          alpha2[k]=next2[k];
        }

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
            emit1[0] = probs1[obs];
            emit2[0] = probs2[obs];
          } else if (obs == obs_block) {
            emit1[1] = 1.0;
            emit2[1] = 1.0;
          } else if (obs == obs_block + 1) {
            emit1[2] = 1.0;
            emit2[2] = 1.0;
          } else if (obs >= obs_offset_6 && obs < obs_offset_6 + obs_block) {
            int idx = obs - obs_offset_6;
            emit1[6] = probs1[idx];
            emit2[6] = probs2[idx];
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
          alpha1[k_state] *= emit1[k_state];
          alpha2[k_state] *= emit2[k_state];
        }

      }

      double L1=0.0;
      double L2=0.0;

      for(int k=0;k<S;k++){
        L1+=alpha1[k];
        L2+=alpha2[k];
      }

      double Li =
        mixing_prop[0]*L1 +
        mixing_prop[1]*L2;

      logL[j] = counts[j]*std::log(Li);

    }

  }

};


// [[Rcpp::export]]
double loglikelihood_ultrafast(
    const NumericMatrix& y,
    const NumericVector& first,
    const NumericMatrix& combins,
    const NumericMatrix& s,
    const NumericMatrix& e,
    const NumericMatrix& i,
    const NumericMatrix& lambda,
    const NumericMatrix& epsilon1,
    const NumericMatrix& epsilon2,
    const NumericMatrix& epsilon3,
    const NumericVector& mixing_prop,
    const NumericVector& probs1,
    const NumericVector& probs2,
    const NumericVector& counts,
    int L
){

  int N = y.nrow();

  NumericVector logL(N);

  UltraFastWorker worker(
      y,first,combins,s,e,i,lambda,
      epsilon1,epsilon2,epsilon3,
      mixing_prop,probs1,probs2,counts,
      L,
      logL);

  parallelFor(0,N,worker);

  return sum(logL);

}

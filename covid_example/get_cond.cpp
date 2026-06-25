#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_cond_cpp(NumericVector thetastar,
                           NumericMatrix alphas,
                           NumericMatrix betas) {
  int p = thetastar.size();
  int B = alphas.ncol();

  NumericVector res(B);
  NumericVector logx(p);

  for (int i = 0; i < p; ++i)
    logx[i] = std::log(thetastar[i]);

  for (int b = 0; b < B; ++b) {
    double sum = 0.0;
    for (int i = 0; i < p; ++i) {
      double a = alphas(i, b);
      double be = betas(i, b);
      // log Gamma density: a*log(be) - lgamma(a) + (a-1)*logx - be*x
      sum += a * std::log(be) -
             R::lgammafn(a) +
             (a - 1.0) * logx[i] -
             be * thetastar[i];
    }
    res[b] = sum;
  }
  return res;
}
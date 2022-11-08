#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector cppsample(int n){
  Rcpp::IntegerVector I = Rcpp::sample(n, n);
  return I;
}

// [[Rcpp::export]]
Rcpp::NumericVector grad_calc_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y,
                                  double alpha, double beta, double gamma, double rho){
  double N = x.length();
  
  double grad_alpha =  2/N*sum((y-gamma+(gamma-rho)/(1+Rcpp::exp(beta*x-alpha)))*((gamma-rho)*Rcpp::exp(beta*x-alpha))/(Rcpp::pow(Rcpp::exp(beta*x-alpha)+1, 2)));
  
  double grad_beta = -2/N*Rcpp::sum((y-gamma+(gamma-rho)/(1+Rcpp::exp(beta*x-alpha)))*((gamma-rho)*x*Rcpp::exp(beta*x-alpha))/(Rcpp::pow(Rcpp::exp(beta*x-alpha)+1, 2)));
  
  double grad_gamma = 2/N*Rcpp::sum((y-gamma+(gamma-rho)/(1+Rcpp::exp(beta*x-alpha)))*(1/(1+Rcpp::exp(beta*x-alpha))-1));
  
  double grad_rho = -2/N*Rcpp::sum((y-gamma+(gamma-rho)/(1+Rcpp::exp(beta*x-alpha)))*(1/(1+Rcpp::exp(beta*x-alpha))));
  
  Rcpp::NumericVector v(4);
  v[0] = grad_alpha; v[1] = grad_beta; v[2] = grad_gamma; v[3] = grad_rho;
  return v;
}


// [[Rcpp::export]]
Rcpp::NumericVector SG_cpp(Rcpp::NumericVector par,
                           int N,
                           Rcpp::NumericVector x,
                           Rcpp::NumericVector y,
                           double gamma,
                           int maxiter = 50,
                           double epsilon = 1e-6){
  for(int k = 0; k < maxiter; k++){
    Rcpp::IntegerVector samp = cppsample(N);
    for(int j = 0; j < N; ++j){
      Rcpp::NumericVector par0 = par;
      int i = samp[j];
      std::cout << j << std::endl;
      std::cout << i << std::endl;
      std::cout << x[i] << std::endl;
      std::cout << y[i] << std::endl;
      std::cout << par0[0] << " " << par0[1] << " " << par0[2] << " " << par0[3] << std::endl;
      Rcpp::NumericVector dummy = grad_calc_cpp(x[i], y[i], par0[0], par0[1], par0[2], par0[3]);
      std::cout << dummy << std::endl;
      Rcpp::NumericVector par = par0 - gamma* grad_calc_cpp(x[i], y[i], par0[0], par0[1], par0[2], par0[3]);
      if(sum(Rcpp::pow(par-par0, 2))<epsilon) {break;}
    }
  }
  return par;
}

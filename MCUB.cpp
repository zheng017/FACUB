#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <R_ext/Applic.h> 

using namespace arma;
using namespace Rcpp;

mat ensure_mat(SEXP x) {
  if (Rf_isMatrix(x)) return as<mat>(x);
  vec v = as<vec>(x);
  return mat(v.memptr(), v.n_elem, 1);
}

inline double log1pexp(double x) {
  if (x > 35) return x; 
  if (x < -10) return std::exp(x);
  return std::log1p(std::exp(x));
}

inline mat sigmoid(const mat& x) {
  return 1.0 / (1.0 + exp(-x));
}

double RMSE(const vec& a, const vec& b) {
  return std::sqrt(mean(square(a - b)));
}


struct DataForB {
  const mat& X, gamma, R_mat;
  const mat& beta; 
  double m;
  int n, p, q;
};


struct DataForBeta {
  const mat& X, gamma, R_mat;
  const vec& B; // Fixed reference
  double m;
  int n, p, q;
};


double fn_B(int n_par, double *par, void *ex) {
  DataForB *d = (DataForB*)ex;
  vec B_curr(par, n_par, false);
  

  mat ll = repmat(B_curr.t(), d->n, 1) + d->X * d->beta.t();
  
  mat log_term(d->n, d->p);
  for(uword i=0; i<ll.n_elem; ++i) log_term(i) = log1pexp(ll(i));
  
  mat obj = d->gamma % ((d->m - d->R_mat) % ll - (d->m - 1.0) * log_term);
  
  return -accu(obj); 
}

void gr_B(int n_par, double *par, double *gr, void *ex) {
  DataForB *d = (DataForB*)ex;
  vec B_curr(par, n_par, false);
  vec gradient(gr, n_par, false);
  
  mat ll = repmat(B_curr.t(), d->n, 1) + d->X * d->beta.t();
  mat sig_ll = sigmoid(ll);
  

  mat res = d->gamma % ( (d->m - d->R_mat) - (d->m - 1.0) * sig_ll );
  
  gradient = -sum(res, 0).t(); 
}


double fn_Beta(int n_par, double *par, void *ex) {
  DataForBeta *d = (DataForBeta*)ex;
  vec beta_vec(par, n_par, false);
  mat beta_curr = reshape(beta_vec, d->p, d->q);
  
  mat ll = repmat(d->B.t(), d->n, 1) + d->X * beta_curr.t();
  
  mat log_term(d->n, d->p);
  for(uword i=0; i<ll.n_elem; ++i) log_term(i) = log1pexp(ll(i));
  
  mat obj = d->gamma % ((d->m - d->R_mat) % ll - (d->m - 1.0) * log_term);
  
  return -accu(obj);
}

void gr_Beta(int n_par, double *par, double *gr, void *ex) {
  DataForBeta *d = (DataForBeta*)ex;
  vec beta_vec(par, n_par, false);
  mat beta_curr = reshape(beta_vec, d->p, d->q);
  vec gradient(gr, n_par, false);
  
  mat ll = repmat(d->B.t(), d->n, 1) + d->X * beta_curr.t();
  mat sig_ll = sigmoid(ll);
  
  mat common = d->gamma % ( (d->m - d->R_mat) - (d->m - 1.0) * sig_ll );
  
  mat grad_mat = common.t() * d->X;
  
  gradient = -vectorise(grad_mat); 
}


void run_optim(int n_params, vec& params, void* data, optimfn fn, optimgr gr, int maxit) {
  double Fmin = 0;
  int fail = 0;
  int fncount = 0;
  int grcount = 0;
  std::vector<int> mask(n_params, 1);
  
  vmmin(n_params, params.memptr(), &Fmin, fn, gr, maxit, 0, 
        mask.data(), -1e20, 1e-8, 10, data, &fncount, &grcount, &fail);
}

// Main MCUB Function
// [[Rcpp::export]]
List MCUB_cpp(SEXP R_in, SEXP V_in, double m, 
              Nullable<NumericVector> b0_in = R_NilValue,
              Nullable<NumericMatrix> beta0_in = R_NilValue,
              int max_iter = 200, double toler = 1e-4, bool trace = false) {
  
  mat R_mat = ensure_mat(R_in);
  int n = R_mat.n_rows;
  int p = R_mat.n_cols;
  
  mat X;
  int q;
  bool has_covariates = false;
  
  if (Rf_isNull(V_in)) {
    q = 1; 
    X = zeros<mat>(n, 1);
  } else {
    X = ensure_mat(V_in);
    q = X.n_cols;
    has_covariates = true;
  }
  
  // Initialization
  vec pai = ones<vec>(p) * 0.5;
  vec B;
  mat beta;
  
  if (b0_in.isNull()) {
    B = ones<vec>(p);
  } else {
    B = as<vec>(b0_in);
  }
  
  if (beta0_in.isNull()) {
    beta = ones<mat>(p, q);
  } else {
    beta = as<mat>(beta0_in);
  }
  
  // Temporary variables
  vec old_pai = pai;
  vec new_pai = pai;
  vec old_B = B;
  vec new_B = B;
  mat old_beta = beta;
  mat new_beta = beta;
  
  mat gamma(n, p);
  mat lchoose_const(n, p);
  
  for(int i=0; i<n; ++i) {
    for(int j=0; j<p; ++j) {
      lchoose_const(i,j) = R::lchoose(m - 1.0, R_mat(i,j) - 1.0);
    }
  }
  
  int iter = 1;
  
  while (iter <= max_iter) {
    if (trace) Rcout << "Iteration: " << iter << "\n";
    Rcpp::checkUserInterrupt();
    
    // --- E-Step ---
    mat ll = repmat(old_B.t(), n, 1) + X * old_beta.t();
    mat csi = sigmoid(ll);
    
    mat log_pai_mat = repmat(log(old_pai.t() + 1e-10), n, 1);
    mat log_1_pai_mat = repmat(log(1.0 - old_pai.t() + 1e-10), n, 1);
    
    mat term1 = log_pai_mat + lchoose_const + (m - R_mat) % log(csi + 1e-10) + (R_mat - 1.0) % log(1.0 - csi + 1e-10);
    
    mat term2 = log_1_pai_mat - std::log(m);
    
    mat num = exp(term1);
    mat den = num + exp(term2);
    gamma = num / den;
    
    // --- M-Step ---
    
    // 1. Update pai
    new_pai = mean(gamma, 0).t(); 
    
    // 2. Update B 
    DataForB d_B = {X, gamma, R_mat, old_beta, m, n, p, q}; 
    run_optim(p, new_B, &d_B, fn_B, gr_B, 100);
    
    // 3. Update Beta if covariates exist
    if (has_covariates) {
      DataForBeta d_Beta = {X, gamma, R_mat, new_B, m, n, p, q};
      vec beta_vec = vectorise(old_beta);
      run_optim(p*q, beta_vec, &d_Beta, fn_Beta, gr_Beta, 100);
      new_beta = reshape(beta_vec, p, q);
    }
    
    // Convergence Check
    double tol1 = RMSE(new_pai, old_pai);
    double tol2 = RMSE(new_B, old_B);
    double tol3 = norm(new_beta - old_beta, "fro");
    
    if (trace) {
      Rcout << "  Tol Pai: " << tol1 << " B: " << tol2 << " Beta: " << tol3 << "\n";
    }
    
    if (tol1 < toler && tol2 < toler && tol3 < toler) break;
    
    old_pai = new_pai;
    old_B = new_B;
    old_beta = new_beta;
    iter++;
  }
  
  mat ll_final = repmat(new_B.t(), n, 1) + X * new_beta.t();
  
  double tol1 = RMSE(new_pai, old_pai);
  double tol2 = RMSE(new_B, old_B);
  double tol3 = norm(new_beta - old_beta, "fro");
  
  return List::create(
    Named("B") = new_B,
    Named("beta") = new_beta,
    Named("pai") = new_pai,
    Named("iter") = iter - 1,
    Named("tol1") = tol1,
    Named("tol2") = tol2,
    Named("tol3") = tol3,
    Named("csi") = sigmoid(ll_final)
  );
}

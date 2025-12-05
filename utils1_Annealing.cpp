// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <R_ext/Applic.h> // Contains vmmin (BFGS)

using namespace arma;
using namespace Rcpp;

// =============================================================================
// Robust Input Handling
// =============================================================================
mat ensure_mat(SEXP x) {
  if (Rf_isMatrix(x)) return as<mat>(x);
  vec v = as<vec>(x);
  return mat(v.memptr(), v.n_elem, 1);
}

// =============================================================================
// Math Helpers
// =============================================================================
inline double log1pexp(double x) {
  if (x > 35) return x; 
  if (x < -10) return std::exp(x);
  return std::log1p(std::exp(x));
}

inline mat sigmoid(const mat& x) {
  return 1.0 / (1.0 + exp(-x));
}

// =============================================================================
// DATA STRUCTURES FOR VMMIN
// =============================================================================

struct GlobalData {
  const mat& X, mu, sigma, pai, alpha, R_mat;
  double m;
  int n_factors, p, q;
};

struct RowData {
  const rowvec& X_i, mu_i, alpha_i, R_i; 
  const mat& Lambda, beta;
  const rowvec& B;
  const rowvec& sigma_i_fixed; 
  double m;
  int p;
  double beta_kl; // --- MODIFICATION 1: Added beta_kl for Annealing ---
};

// =============================================================================
// VMMIN WRAPPERS
// =============================================================================

// --- 1. Global Optimization (Unchanged, as Lambda/B/beta are not in KL term) ---

double global_fn(int n, double *par, void *ex) {
  GlobalData *d = (GlobalData*)ex;
  vec params(par, n, false); 
  
  mat Lambda = reshape(params.subvec(0, d->p * d->n_factors - 1), d->p, d->n_factors);
  vec B = params.subvec(d->p * d->n_factors, d->p * d->n_factors + d->p - 1);
  mat beta = reshape(params.subvec(d->p * d->n_factors + d->p, n - 1), d->p, d->q);
  
  mat ll = repmat(B.t(), d->mu.n_rows, 1) + d->X * beta.t() + d->mu * Lambda.t();
  mat e_mat = ll + 0.5 * (d->sigma * square(Lambda).t());
  
  mat log_term(e_mat.n_rows, e_mat.n_cols);
  for(uword i=0; i<e_mat.n_elem; ++i) log_term(i) = log1pexp(e_mat(i));
  
  mat term1 = d->alpha % (d->m - d->R_mat) % ll;
  mat term2 = -d->alpha * (d->m - 1.0) % log_term;
  
  return -(accu(term1) + accu(term2)); 
}

void global_gr(int n, double *par, double *gr, void *ex) {
  GlobalData *d = (GlobalData*)ex;
  vec params(par, n, false);
  vec gradient(gr, n, false); 
  
  mat Lambda = reshape(params.subvec(0, d->p * d->n_factors - 1), d->p, d->n_factors);
  vec B = params.subvec(d->p * d->n_factors, d->p * d->n_factors + d->p - 1);
  mat beta = reshape(params.subvec(d->p * d->n_factors + d->p, n - 1), d->p, d->q);
  
  mat ll = repmat(B.t(), d->mu.n_rows, 1) + d->X * beta.t() + d->mu * Lambda.t();
  mat e_mat = ll + 0.5 * (d->sigma * square(Lambda).t());
  mat sig_e = sigmoid(e_mat);
  
  mat Common = d->alpha % ((d->m - d->R_mat) - (d->m - 1.0) * sig_e);
  mat Weight = -d->alpha * (d->m - 1.0) % sig_e;
  
  mat grad_Lambda_mat = Common.t() * d->mu + (Weight.t() * d->sigma) % Lambda;
  vec grad_B_vec = sum(Common, 0).t();
  mat grad_beta_mat = Common.t() * d->X;
  
  gradient.subvec(0, d->p * d->n_factors - 1) = vectorise(-grad_Lambda_mat);
  gradient.subvec(d->p * d->n_factors, d->p * d->n_factors + d->p - 1) = -grad_B_vec;
  gradient.subvec(d->p * d->n_factors + d->p, n - 1) = vectorise(-grad_beta_mat);
}

// --- 2. Sigma Row Wrappers (Modified for Annealing) ---

double sigma_fn(int n, double *par, void *ex) {
  RowData *d = (RowData*)ex;
  vec log_sigma(par, n, false);
  rowvec sigma_i = exp(log_sigma.t());
  
  rowvec ll = d->B + d->X_i * d->beta.t() + d->mu_i * d->Lambda.t();
  rowvec e_mat = ll + 0.5 * (sigma_i * square(d->Lambda).t());
  
  // --- MODIFICATION 1: Apply beta_kl ---
  double term1 = -0.5 * d->beta_kl * (accu(sigma_i) - accu(log_sigma));
  
  rowvec log_term(d->p);
  for(int j=0; j<d->p; ++j) log_term(j) = log1pexp(e_mat(j));
  double term2 = -accu( (d->alpha_i * (d->m-1.0)) % log_term );
  
  return -(term1 + term2);
}

void sigma_gr(int n, double *par, double *gr, void *ex) {
  RowData *d = (RowData*)ex;
  vec log_sigma(par, n, false);
  vec gradient(gr, n, false);
  
  rowvec sigma_i = exp(log_sigma.t());
  rowvec ll = d->B + d->X_i * d->beta.t() + d->mu_i * d->Lambda.t();
  rowvec e_mat = ll + 0.5 * (sigma_i * square(d->Lambda).t());
  mat sig_e = sigmoid(e_mat);
  
  rowvec A = -(d->m - 1.0) * d->alpha_i % sig_e;
  
  // --- MODIFICATION 1: Apply beta_kl to prior gradient ---
  // Derivative of -0.5*beta*(sigma - log(sigma)) wrt sigma is -0.5*beta*(1 - 1/sigma)
  rowvec d_vlb_d_sigma = -0.5 * d->beta_kl * (1.0 - 1.0/sigma_i) + A * square(d->Lambda);
  
  // Chain rule for log_sigma parametrization
  rowvec d_vlb_d_log = d_vlb_d_sigma % sigma_i;
  
  gradient = -d_vlb_d_log.t();
}

// --- 3. Mu Row Wrappers (Modified for Annealing) ---

double mu_fn(int n, double *par, void *ex) {
  RowData *d = (RowData*)ex;
  vec mu_i(par, n, false);
  
  rowvec mu_row = mu_i.t();
  rowvec ll = d->B + d->X_i * d->beta.t() + mu_row * d->Lambda.t();
  rowvec e_mat = ll + 0.5 * (d->sigma_i_fixed * square(d->Lambda).t());
  
  // --- MODIFICATION 1: Apply beta_kl ---
  double term1 = -0.5 * d->beta_kl * accu(square(mu_row));
  
  double term2 = accu( d->alpha_i % ((d->m - d->R_i) % (mu_row * d->Lambda.t())) );
  
  rowvec log_term(d->p);
  for(int j=0; j<d->p; ++j) log_term(j) = log1pexp(e_mat(j));
  double term3 = -accu( (d->alpha_i * (d->m-1.0)) % log_term );
  
  return -(term1 + term2 + term3);
}

void mu_gr(int n, double *par, double *gr, void *ex) {
  RowData *d = (RowData*)ex;
  vec mu_i(par, n, false);
  vec gradient(gr, n, false);
  
  rowvec mu_row = mu_i.t();
  rowvec ll = d->B + d->X_i * d->beta.t() + mu_row * d->Lambda.t();
  rowvec e_mat = ll + 0.5 * (d->sigma_i_fixed * square(d->Lambda).t());
  
  rowvec sum1 = d->alpha_i % ( (d->m - d->R_i) - (d->m - 1.0) * sigmoid(e_mat) );
  
  // --- MODIFICATION 1: Apply beta_kl to prior gradient ---
  // Derivative of -0.5*beta*mu^2 is -beta*mu
  rowvec grad = -d->beta_kl * mu_row + sum1 * d->Lambda;
  
  gradient = -grad.t();
}

// =============================================================================
// Helper to Call vmmin
// =============================================================================
void run_optim(int n_params, vec& params, void* data, optimfn fn, optimgr gr, int maxit) {
  double Fmin = 0;
  int fail = 0;
  int fncount = 0;
  int grcount = 0;
  
  std::vector<int> mask(n_params, 1);
  
  vmmin(n_params, params.memptr(), &Fmin,
        fn, gr,
        maxit, 0, 
        mask.data(), 
        -1e20, 1e-6, 
        10, 
        data, &fncount, &grcount, &fail);
}

// =============================================================================
// Main Function
// =============================================================================
// [[Rcpp::export]]
List FAVA_KLAnnealing_cpp(SEXP R_in, SEXP V_in, double m, int n_factors, 
                          List initials, int maxit = 200, bool trace = false) {
  
  mat R_mat = ensure_mat(R_in);
  int n = R_mat.n_rows;
  int p = R_mat.n_cols;
  
  mat X;
  int q;
  if (Rf_isNull(V_in)) { q = 1; X = zeros<mat>(n, 1); } 
  else { X = ensure_mat(V_in); q = X.n_cols; }
  
  mat mu = ensure_mat(initials["mu"]);
  mat sigma = ensure_mat(initials["sigma"]);
  mat pai = ensure_mat(initials["pai"]);
  mat alpha = ensure_mat(initials["alpha"]);
  vec B = as<vec>(initials["B"]);
  mat beta = ensure_mat(initials["beta"]);
  mat Lambda = ensure_mat(initials["Lambda"]);
  
  mat new_mu = mu;
  mat new_sigma = sigma;
  mat new_pai = pai;
  mat new_alpha = alpha;
  vec new_B = B;
  mat new_beta = beta;
  mat new_Lambda = Lambda;
  
  double cur_VLB = -1e6; 
  double eps = 1e-4;
  int iter = 1;
  double diff = 1e10; 
  
  mat log_choose_const(n, p);
  for(int i=0; i<n; ++i) 
    for(int j=0; j<p; ++j) 
      log_choose_const(i,j) = R::lchoose(m - 1.0, R_mat(i,j) - 1.0);
  
  // --- MODIFICATION 1: Annealing Parameters ---
  int anneal_steps = 50; 
  
  while ( (diff > eps * (std::abs(cur_VLB) + eps)) && (iter <= 100) ) {
    
    // --- MODIFICATION 1: Calculate beta_kl ---
    double beta_kl = (iter <= anneal_steps) ? (double)iter/anneal_steps : 1.0;
    
    if (trace) Rcout << "Iteration: " << iter << " (KL Weight: " << beta_kl << ")\n";
    Rcpp::checkUserInterrupt();
    
    // --- MODIFICATION 2: MAP Estimation for Pai (Smoothing) ---
    // Using Beta(2,2) prior -> (sum + 1) / (n + 2)
    rowvec smooth_pai = (sum(new_alpha, 0) + 1.0) / (n + 2.0);
    new_pai = repmat(smooth_pai, n, 1);
    
    // if (beta_kl >= 0.95) {
    //   rowvec smooth_pai = (sum(new_alpha, 0) + 1.0) / (n + 2.0);
    //   smooth_pai = clamp(smooth_pai, 0.001, 0.99); 
    //   new_pai = repmat(smooth_pai, n, 1);
    // } 
    
    
    // --- 1. Global Optimization (using vmmin) ---
    vec params = join_vert(vectorise(new_Lambda), new_B);
    params = join_vert(params, vectorise(new_beta));
    
    GlobalData g_data = {X, new_mu, new_sigma, new_pai, new_alpha, R_mat, m, n_factors, p, q};
    run_optim(params.n_elem, params, &g_data, global_fn, global_gr, maxit);
    
    new_Lambda = reshape(params.subvec(0, p * n_factors - 1), p, n_factors);
    new_B = params.subvec(p * n_factors, p * n_factors + p - 1);
    new_beta = reshape(params.subvec(p * n_factors + p, params.n_elem - 1), p, q);
    
    // --- Inner Loop ---
    double delta_alpha_req = 1e-3;
    int p_iter = 1;
    bool inner_converged = false;
    
    while (!inner_converged && p_iter < 10) {
      
      // Save previous values for robust stopping
      mat prev_alpha = new_alpha;
      mat prev_mu = new_mu;
      mat prev_sigma = new_sigma;
      
      // 3a. Update Alpha
      mat ll = repmat(new_B.t(), n, 1) + X * new_beta.t() + new_mu * new_Lambda.t();
      mat e_mat = ll + 0.5 * (new_sigma * square(new_Lambda).t());
      mat log_term_e(n, p);
      for(uword k=0; k<e_mat.n_elem; ++k) log_term_e(k) = log1pexp(e_mat(k));
      mat term_logit_pai(n, p);
      for(uword k=0; k<new_pai.n_elem; ++k) {
        double val = std::max(1e-8, std::min(new_pai(k), 1.0 - 1e-8));
        term_logit_pai(k) = std::log(val / (1.0 - val));
      }
      mat log_odds = term_logit_pai + log_choose_const + (m - R_mat) % ll - (m - 1.0) * log_term_e + std::log(m);
      new_alpha = sigmoid(log_odds);
      
      // 3b. Update Sigma
      for(int i=0; i<n; ++i) {
        vec log_sigma_i = log(new_sigma.row(i).t());
        RowData s_data = {X.row(i), new_mu.row(i), new_alpha.row(i), R_mat.row(i), 
                          new_Lambda, new_beta, new_B.t(), new_sigma.row(i), m, p, beta_kl}; // Pass beta_kl
        
        run_optim(n_factors, log_sigma_i, &s_data, sigma_fn, sigma_gr, 20);
        new_sigma.row(i) = exp(log_sigma_i).t();
      }
      
      // 3c. Update Mu
      for(int i=0; i<n; ++i) {
        vec mu_i = new_mu.row(i).t();
        RowData m_data = {X.row(i), new_mu.row(i), new_alpha.row(i), R_mat.row(i), 
                          new_Lambda, new_beta, new_B.t(), new_sigma.row(i), m, p, beta_kl}; // Pass beta_kl
        
        run_optim(n_factors, mu_i, &m_data, mu_fn, mu_gr, 20);
        new_mu.row(i) = mu_i.t();
      }
      
      // --- MODIFICATION 3: Robust Inner Stopping ---
      double max_d_alpha = max(vectorise(abs(new_alpha - prev_alpha)));
      double max_d_mu    = max(vectorise(abs(new_mu - prev_mu)));
      double max_d_sigma = max(vectorise(abs(new_sigma - prev_sigma)));
      
      // Check if ALL parameters stabilized
      if (max_d_alpha < delta_alpha_req && max_d_mu < 1e-3 && max_d_sigma < 1e-3) {
        inner_converged = true;
      }
      
      p_iter++;
    }
    
    // 4. Convergence Check (VLB)
    mat ll = repmat(new_B.t(), n, 1) + X * new_beta.t() + new_mu * new_Lambda.t();
    mat e_mat = ll + 0.5 * (new_sigma * square(new_Lambda).t());
    
    // --- MODIFICATION 1: Apply beta_kl to VLB check for consistency ---
    double fun1 = -0.5 * beta_kl * (accu(new_sigma) + accu(square(new_mu)) - accu(log(new_sigma)));
    
    mat eps_mat = ones<mat>(n, p) * 1e-8;
    mat t2_a = (new_alpha + eps_mat) % log((new_pai + eps_mat) / (new_alpha + eps_mat));
    mat t2_b = abs(1.0 - new_alpha - eps_mat) % log((1.0 - new_pai - eps_mat) / abs(1.0 - new_alpha - eps_mat));
    double fun2 = accu(t2_a + t2_b);
    mat log_term_e(n, p);
    for(uword k=0; k<e_mat.n_elem; ++k) log_term_e(k) = log1pexp(e_mat(k));
    double fun3 = accu( new_alpha % ( log_choose_const + (m - R_mat) % ll - (m - 1.0) * log_term_e ) ) 
      - accu((1.0 - new_alpha) * std::log(m));
    
    double new_VLB = fun1 + fun2 + fun3;
    diff = std::abs(new_VLB - cur_VLB);
    cur_VLB = new_VLB;
    iter++;
  }
  
  mat ll = repmat(new_B.t(), n, 1) + X * new_beta.t() + new_mu * new_Lambda.t();
  mat e_mat = ll + 0.5 * (new_sigma * square(new_Lambda).t());
  mat log_term_ll(n, p); for(uword k=0; k<ll.n_elem; ++k) log_term_ll(k) = log1pexp(ll(k));
  mat log_term_e(n, p);  for(uword k=0; k<e_mat.n_elem; ++k) log_term_e(k) = log1pexp(e_mat(k));
  double cur_log = accu( new_alpha % ( log_choose_const + (m - R_mat)%ll - (m-1.0)*log_term_ll ) - (1.0 - new_alpha) * std::log(m) );
  double EE = accu( new_alpha % ( log_choose_const + (m - R_mat)%ll - (m-1.0)*log_term_e ) - (1.0 - new_alpha) * std::log(m) );
  
  return List::create(
    Named("VLB") = cur_VLB, Named("EE") = EE, Named("lob") = cur_log, Named("EN2") = 2.0 * cur_log - 2.0 * EE,
          Named("iter") = iter - 1,
          Named("lvs") = List::create(Named("alpha") = new_alpha, Named("mu") = new_mu, Named("sigma") = new_sigma),
                Named("params") = List::create(Named("pai") = new_pai.row(0), Named("B") = new_B, Named("beta") = new_beta, Named("Lambda") = new_Lambda),
                Named("ll") = sigmoid(ll), Named("e.mat") = sigmoid(e_mat)
  );
}

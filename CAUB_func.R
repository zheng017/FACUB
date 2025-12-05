library(mvnfast)

MSE = function(a, b) {
  mean((a-b)^2)
}

RMSE = function(a, b) {
  sqrt(MSE(a, b))
}

prob_u = function(m, alpha) {
  prob = numeric(m)
  for (r in 1:m) {
    a = (r-1)/m
    b = r/m
    prob[r] = pbeta(b, alpha, alpha) - pbeta(a, alpha, alpha)
  }
  prob
}

generate_data_CAUB = function(n, m, pai, csi, alpha) {
  # generate n samples of a ordinal variable, distributed from CAUB
  # pai, csi, alpha are scalars, rather than vectors
  # the output is a n dimensional vector of samples
  z = sample(0:1, n, replace = TRUE, prob = c(1-pai, pai))
  z_idx = which(z==1)
  res = numeric(n)
  res[z_idx] = rbinom(length(z_idx), m-1, 1-csi) + 1
  res[setdiff(1:n, z_idx)] = sample(1:m, n-length(z_idx), replace = TRUE, prob = prob_u(m, alpha))
  return(res)
}

generate_data_FACAUB = function(n, p, q, d, m, b0, pai_param, zeta, type1 = "uniform", type2 = "identity") {
  a = rep(0, n)
  b = rep(b0, p)
  # generate loadings matrix: Lambda  pxd
  if (type1 == "uniform") {
    Lambda = matrix(runif(p*d, -2, 2), p, d)
  } else if (type1 == "gaussian") {
    Lambda = matrix(rnorm(p*d), p, d)
  }
  # generate factors: f nxd
  f = matrix(rnorm(n*d), n, d)
  if (q == 0) {
    x1 = NULL
    beta = NULL
    csi = plogis( matrix(a, n, p) + matrix(b, n, p, byrow = TRUE) + f %*% t(Lambda) )
  } else {
    beta = matrix(rnorm(p*q), p, q)
    if (type2 == "identity") {
      x = matrix(rnorm(n*q), n, q)
    } else if (type2 == "exchange") {
      rho = 0.3
      Sigma = matrix(rho, q, q)
      diag(Sigma) = rep(1, q)
      x = rmvn(n, rep(0, q), Sigma)
    }
    csi = plogis( matrix(a, n, p) + matrix(b, n, p, byrow = TRUE) + x1 %*% beta + f %*% t(Lambda) )
  }
  
  pai = runif(p, pai_param[1], pai_param[2])
  z = R = matrix(NA, n, p)
  for (j in 1:p) {
    z[, j] = sample(0:1, n, replace = TRUE, prob = c(1-pai[j], pai[j]))
    z_idx = which(z[, j] == 1)
    R[z_idx, j] = rbinom(length(z_idx), m-1, prob = 1-csi[z_idx, j]) + 1
    z0_idx = setdiff(1:n, z_idx)
    R[z0_idx, j] = sample(1:m, length(z0_idx), replace = TRUE, prob = prob_u(m, zeta[j]))
  }
  pai_mat = matrix(pai, n, p, byrow = TRUE)
  prob = pai_mat*prob_feeling(m, R, csi) + (1-pai_mat)*prob_uncertainty(m, R, matrix(zeta, n, p, byrow = TRUE))
  return(list(R = R,
              A = a,
              B = b,
              Lambda = Lambda,
              f = f,
              csi = csi,
              pai = pai,
              zeta = zeta,
              prob = prob))
}



generate_data_FACAUB2 = function(n, p, q, d, m, b0, pai_param, zeta, type1 = "uniform", type2 = "identity") {
  a = rep(0, n)
  b = rep(b0, p)
  # generate loadings matrix: Lambda  pxd
  Lambda = matrix(NA, p, d)
  if (type1 == "uniform") {
    for (i in 1:p) {
      Lambda[i, ] = runif(d, -2, 2)
    }
  } else if (type1 == "gaussian") {
    Lambda[i, ] = rnorm(d, 0, 1)
  }
  # generate factors: f nxd
  f = matrix(NA, n, d)
  for (i in 1:n) {
    f[i, ] = rnorm(d, 0, 1)
  }
  
  if (q == 0) {
    x1 = NULL
    beta = NULL
    csi = plogis( matrix(a, n, p) + matrix(b, n, p, byrow = TRUE) + f %*% t(Lambda) )
  } else {
    beta = matrix(rnorm(p*q), p, q)
    if (type2 == "identity") {
      x = matrix(rnorm(n*q), n, q)
    } else if (type2 == "exchange") {
      rho = 0.3
      Sigma = matrix(rho, q, q)
      diag(Sigma) = rep(1, q)
      x = rmvn(n, rep(0, q), Sigma)
    }
    csi = plogis( matrix(a, n, p) + matrix(b, n, p, byrow = TRUE) + x1 %*% beta + f %*% t(Lambda) )
  }
  
  pai = runif(p, pai_param[1], pai_param[2])
  z = R = matrix(NA, n, p)
  nk = numeric(n)
  for (i in 1:n) {
    z[i, ] = rbinom(p, size = 1, prob = pai)
    idx0 = which(z[i, ] == 0)
    nk[i] = length(idx0)
    for (j in 1:p) {
      R[i, j] = rbinom(1, size=m-1, prob=1-csi[i, j]) + 1
    }
    if (nk[i] != 0) {
      for (k in 1:nk[i]) {
        R[i, idx0[k]] = sample(1:m, 1, prob = prob_u(m, zeta[idx0[k]]))
      }
    }
  }
  
  pai_mat = matrix(pai, n, p, byrow = TRUE)
  prob = pai_mat*prob_feeling(m, R, csi) + (1-pai_mat)*prob_uncertainty(m, R, matrix(zeta, n, p, byrow = TRUE))
  return(list(R = R,
              A = a,
              B = b,
              Lambda = Lambda,
              f = f,
              csi = csi,
              pai = pai,
              zeta = zeta,
              prob = prob))
}

generate_data_CAUB_covariates = function(n, m, p, q, pai) {
  # generate n samples from a ordinal variable, distributed from CAUB, 
  # with csi and alpha modeled in a regression model with
  # logit(csi) = x1 %*% beta1,
  # log(alpha) = x2 %*% beta2.
  # pai is a scalar, as we only consider one ordinal variable, and do not consider pai_i,
  # that is, individual-specific pai.
  x1 = cbind(rep(1, n), rmvn(n, rep(0, p-1), diag(p-1)))
  x2 = cbind(rep(1, n), rmvn(n, rep(0, q-1), diag(q-1)))
  beta1 = runif(p, 0, 2)
  beta2 = runif(q, -.5, .5)
  csi = plogis(x1 %*% beta1)
  alpha = exp(x2 %*% beta2)
  
  z = sample(0:1, n, replace = TRUE, prob = c(1-pai, pai))
  z_idx = which(z==1)
  res = numeric(n)
  res[z_idx] = rbinom(length(z_idx), m-1, 1-csi[z_idx]) + 1
  z0_idx = setdiff(1:n, z_idx)
  for (i in z0_idx) {
    res[i] = sample(1:m, 1, replace = TRUE, prob = prob_u(m, alpha[i]))
  }
  prob = pai*prob_feeling(m, res, csi) + (1-pai)*prob_uncertainty(m, res, alpha)
  return(list(data = res,
              x1 = x1,
              x2 = x2,
              beta1 = beta1,
              beta2 = beta2,
              csi = csi,
              alpha = alpha,
              pai = pai,
              prob = prob))
}


prob_feeling = function(m, r, csi) {
  # csi can be a vector
  choose(m-1, r-1) * csi^(m-r) * (1-csi)^(r-1)
}

prob_uncertainty = function(m, r, alpha, islog = FALSE) {
  # alpha can be a vector
  out = pbeta(r/m, alpha, alpha) - pbeta((r-1)/m, alpha, alpha)
  out[out<=0] = 1e-08
  if (islog) {
    return(log(out))
  } else {
    return(out)
  }
}


CAUB = function(r, m, iter = 1000, tol = 1e-3) {
  pai_hat = pai_old = 0.5
  csi_hat = csi_old = 0.5
  alpha_hat = alpha_old = 1
  count = 1
  for (i in 1:iter) {
    tau = pai_old * prob_feeling(m, r, csi_old) / (pai_old * prob_feeling(m, r, csi_old) + (1-pai_old) * prob_uncertainty(m, r, alpha_old))
    pai_hat = mean(tau)
    csi_hat = sum(tau*(m-r)) / ((m-1)*sum(tau))
    fn = function(x) {
      -sum((1-tau)*log(pbeta(r/m, x, x)-pbeta((r-1)/m, x, x)))
    }
    alpha_hat = optim(alpha_old, fn, method = "BFGS", lower = 1e-6)$par
    if (abs(pai_hat-pai_old) < tol & abs(csi_hat-csi_old) < tol & abs(alpha_hat - alpha_old) < tol) break
    
    pai_old = pai_hat
    alpha_old = alpha_hat
    csi_old = csi_hat
    count = count + 1
  }
  return(list(pai = pai_hat,
              alpha = alpha_hat,
              csi = csi_hat,
              iter = count))
}

CAUB_covariates = function(r, m, x1, x2, iter = 1000, tol = 1e-3, grad = FALSE, type = "BFGS") {
  p = ncol(x1)
  q = ncol(x2)
  pai_hat = pai_old = 0.5
  beta1_old = rep(1, p)
  beta2_old = runif(q, -.5, .5)
  csi_old = plogis(x1 %*% beta1_old)
  alpha_old = exp(x2 %*% beta2_old)
  
  count = 1
  for (i in 1:iter) {
    tau = pai_old * prob_feeling(m, r, csi_old) / (pai_old * prob_feeling(m, r, csi_old) + (1-pai_old) * prob_uncertainty(m, r, alpha_old))
    pai_hat = mean(tau)
    fn_beta = function(x) {
      -sum( tau * ( (m-r) * (x1 %*% x) - (m-1)*log1p(exp(x1 %*% x)) ) )
    }
    gn_beta = function(x) {
      -colSums(as.vector( tau * ( (m-r) - (m-1)*plogis(x1 %*% x) ) ) * x1)
    }
    beta1_hat = optim(beta1_old, fn_beta, gr = gn_beta, method = "BFGS")$par
    csi_hat = plogis(x1 %*% beta1_hat)
    
    fn_beta2 = function(x) {
      alpha_tmp = exp(x2 %*% x)
      tmp = prob_uncertainty(m, r, alpha_tmp, islog = TRUE)
      -sum( (1-tau) *  tmp )
    }
    # history = list()
    # fn_beta2_debug = function(x) {
    #   val = fn_beta2(x)
    #   history <<- append(history, list(par = x, value = val))
    #   val
    # }
    gn_beta2 = function(x) {
      alpha_tmp = exp(x2 %*% x)
      integral = integral2 = numeric(n)
      for (i in 1:n) {
        integral[i] = integrate(f = function(t) (t^(alpha_tmp[i]-1))*(1-t)^(alpha_tmp[i]-1)*log(t*(1-t)), lower = (r[i]-1)/m, upper = r[i]/m)$value
        integral2[i] = integrate(f = function(t) (t^(alpha_tmp[i]-1))*(1-t)^(alpha_tmp[i]-1)*log(t*(1-t)), lower = 0, upper = 1)$value
      }
      tmp = prob_uncertainty(m, r, alpha_tmp)
      tmp = tmp * beta(alpha_tmp, alpha_tmp)
      tmp2 = beta(alpha_tmp, alpha_tmp)
      -colSums(as.vector( (1-tau) * (integral / tmp - integral2 / tmp2) * alpha_tmp ) * x2)
    }
    if (grad) {
      if (type == "BFGS") {
        beta2_hat = optim(beta2_old, fn_beta2, gr = gn_beta2, method = "BFGS")$par
      } else {
        beta2_hat = optim(beta2_old, fn_beta2, gr = gn_beta2)$par
      }
    } else {
      if (type == "BFGS") {
        beta2_hat = optim(beta2_old, fn_beta2, method = "BFGS")$par
      } else {
        beta2_hat = optim(beta2_old, fn_beta2)$par
      }
      #print(history)
      #print(sapply(history, function(x) is.nan(x$value) || any(is.nan(x$par))))
    }
    alpha_hat = exp(x2 %*% beta2_hat)
    if (abs(pai_hat-pai_old) < tol & norm(csi_hat - csi_old, "2") < tol & norm(alpha_hat - alpha_old, "2") < tol) break
    
    pai_old = pai_hat
    alpha_old = alpha_hat
    csi_old = csi_hat
    count = count + 1
  }
  prob = pai_hat*prob_feeling(m, r, csi_hat) + (1-pai)*prob_uncertainty(m, r, alpha_hat)
  return(list(pai = pai_hat,
              beta1 = beta1_hat,
              beta2 = beta2_hat,
              alpha = alpha_hat,
              csi = csi_hat,
              prob = prob,
              iter = count))
}

MCUB = function(R, V=NULL, m, b0=rep(1,ncol(R)), beta0=matrix(1,ncol(R), ifelse(is.null(V), 1, ncol(V))), max.iter=200, toler=1e-4, trace=FALSE) {
  # 使用EM算法来估计模型参数
  n = nrow(R)
  p = ncol(R)
  if (is.null(V)) {
    q = 1
    X = matrix(0, n, q)
  } else {
    X = V
    q = ncol(X)
  }
  
  # Initialization
  pai = new.pai = old.pai = rep(0.5, p)
  B = new.B = old.B = b0
  beta = new.beta = old.beta = beta0
  
  
  Q_function = function(R, gamma, B, beta, pai) {
    # E-step 完全数据的log似然求条件期望
    # gamma为E[z|data,previous]
    # 此函数也可视为关于优化B, beta的目标函数
    # 这里输入的beta是pq*1维的向量
    beta = matrix(beta, p, q)  # 把pq*1维向量转化成p*q矩阵
    ll = matrix(B, n, p, byrow = TRUE) + X %*% t(beta)
    #res = julia_call("compute_res_jl", n, m, R, pai, ll, gamma)
    res = gamma*(log(matrix(pai, n, p, byrow = TRUE)) + log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(ll))) +
      (1-gamma)*(log(matrix(1-pai, n, p, byrow = TRUE)) - log(m))
    return(sum(res))
  }
  
  gr_B = function(R, gamma, B, beta, pai) {
    beta = matrix(beta, p, q)
    ll = matrix(B, n, p, byrow = TRUE) + X %*% t(beta)
    res = gamma*( m-R - (m-1)*exp(ll)/(1+exp(ll)) )
    return(colSums(res))
  }
  
  gr_beta = function(R, gamma, B, beta, pai) {
    beta = matrix(beta, p, q)
    ll = matrix(B, n, p, byrow = TRUE) + X %*% t(beta)
    b2 = NULL
    for (w in 1:q) {
      b2 = c(b2, gamma*( sweep(m-R - (m-1)*exp(ll)/(1+exp(ll)), 1, X[, w], "*") ))
    }
    b2 = matrix(b2, n, p*q)
    return(colSums(b2))
  }
  
  # EM algorithm
  iter = 1; b.cur.logfunc = beta.cur.logfunc = -1e6
  while (iter <= max.iter) {
    # E-step
    # 计算gamma=E[z|data,previous]
    ll = matrix(old.B, n, p, byrow = TRUE) + X %*% t(old.beta)
    csi = plogis(ll)
    fenzi = log(matrix(old.pai, n, p, byrow = TRUE)) + log(choose(m-1, R-1)) + (m-R)*log(csi) + (R-1)*log(1-csi+1e-8)
    fenmu = exp(fenzi) + (matrix(1-old.pai, n, p, byrow = TRUE)) / m
    gamma = exp(fenzi) / fenmu
    
    # M-step
    # 更新pai
    new.pai = colMeans(gamma)
    # 更新B
    qq = try(optim(old.B, fn = Q_function, gr = gr_B, R = R, gamma = gamma, beta = c(new.beta), pai = new.pai, method = "BFGS", control = list(trace=0, fnscale=-1, maxit=100)), silent = TRUE)
    if ("try-error" %in% class(qq)) {
      new.B = old.B
    } else {
      if (iter > 1 && b.cur.logfunc > qq$value) {
        if (trace) {cat("Optimization of model coefs B did not improve on iteration step ", iter, "\n")}
        new.B = old.B
      } else {
        if (trace) {cat("Model parameters updated", "\n")}
        new.B = qq$par
        if (qq$convergence != 0) {
          if (trace) {cat("Optimization of model coefs B did not converge on iteration step ", iter, "\n")}
        }
      }
    }
    
    if (!is.null(V)) {
      # 更新beta
      qq = try(optim(c(old.beta), fn = Q_function, gr = gr_beta, R = R, gamma = gamma, B = new.B, pai = new.pai, method = "BFGS", control = list(trace=0, fnscale=-1, maxit=100)), silent = TRUE)
      if ("try-error" %in% class(qq)) {
        new.beta = old.beta
      } else {
        if (iter > 1 && beta.cur.logfunc > qq$value) {
          if (trace) {cat("Optimization of model coefs beta did not improve on iteration step ", iter, "\n")}
          new.beta = old.beta
        } else {
          if (trace) {cat("Model parameters updated", "\n")}
          new.beta = qq$par
          new.beta = matrix(new.beta, p, q)
          if (qq$convergence != 0) {
            if (trace) {cat("Optimization of model coefs beta did not converge on iteration step ", iter, "\n")}
          }
        }
      }
    }
    # 终止准则
    tol1 = RMSE(new.pai, old.pai)
    tol2 = RMSE(new.B, old.B)
    tol3 = norm(new.beta - old.beta, "F")
    if ( (tol1 < toler) && (tol2 < toler) && (tol3 < toler) ) break;
    
    old.pai = new.pai
    old.B = new.B
    old.beta = new.beta
    b.cur.logfunc = beta.cur.logfunc = Q_function(R, gamma, new.B, c(new.beta), new.pai)
    iter = iter + 1
  }
  
  ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta)
  
  return(list(B = new.B, beta = new.beta, pai = new.pai, iter = iter-1, tol1 = tol1, tol2 = tol2, tol3 = tol3, csi = plogis(ll)))
}


check_R = function(R, m) {
  if (m %% 2 == 1) {
    left = 1:((m-1)/2)
    right = ((m+1)/2+1):m
  } else if (m %% 2 == 0) {
    left = 1:(m/2)
    right = (m/2+1):m
  }
  reorder = FALSE
  freq = table(R)
  if (sum(freq[left]) >= sum(freq[right])) {
    R = m+1-R
    reorder = TRUE
  }
  return(list(R = R,
              reorder = reorder))
}

# no covariate
FACAUB = function(R, x1=NULL, x2=NULL, m, n_factors, maxit = 200, trace = FALSE, init = NULL) {
  n = nrow(R)
  p = ncol(R)
  # if covariates about csi included
  if (is.null(x1)) {
    q = 1
    x1 = matrix(0, n, q)
  } else {
    q = ncol(x1)
  }
  
  # initialization of parameters and variational parameters
  R_scale = scale(log2(1+R), center = TRUE, scale = TRUE)
  re = svd(R_scale, n_factors, n_factors)
  if (n_factors == 1) {
    mu = new_mu = re$u * (re$d[1])
  } else {
    mu = new_mu = re$u %*% diag(re$d[1:n_factors])  # nxn_factors dimension
  }
  Lambda = new_Lambda = re$v
  B = new_B = rep(1, p)
  beta = new_beta = matrix(1, p, q)
  sigma = new_sigma = matrix(0.01, n, n_factors)  # nxn_factors dimension
  pai = new_pai = rep(0.5, p)
  alpha = new_alpha = matrix(pai, n, p, byrow = TRUE)
  if (!is.null(init)) {
    mu=new.mu=init$mu
    sigma=new.sigma=init$sigma
    B=new.B=init$B
    beta=new.beta=init$beta
    Lambda=new.Lambda=init$Lambda
    pai=new.pai=init$pai
    alpha=new.alpha=init$alpha
  }
  # parameters of discrete beta distribution in adjusted uncertainty term
  # first assume that no covariates included in zeta, just as parameters pai.
  zeta = new_zeta = rep(1, p)
  
  # VA iteration about variational lower bound
  cur.VLB = -1e6; iter = 1; ratio = 10; diff = 1e5; eps = 1e-4; max.iter = 100;
  m.cur.logfunc = -1e6; b.cur.logfunc = -1e6; d.cur.logfunc = -1e6; s.cur.logfunc = -1e6;
  
  while ( (diff > eps*(abs(cur.VLB)+eps)) && iter <= max.iter ) {
    if (trace) cat("Iteration: ", iter, "\n")
    
    # variational lower bound
    VLB = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
      new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
      new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
      new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
      model_coefs = model_coefs[-(1:(p*n_factors))]
      new_B = model_coefs[1:p]
      model_coefs = model_coefs[-(1:p)]
      new_beta = matrix(model_coefs[1:(p*q)], p, q)
      
      new_zeta = zeta
      new_pai = pai
      new_alpha = alpha
      
      ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
      e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
      
      fun1 = -0.5*( sum(new_sigma) + sum(new_mu^2) - sum(log(new_sigma)) )
      fun2 = (new_alpha+1e-8) * log( matrix(new_pai+1e-8, n, p, byrow = TRUE) / (new_alpha+1e-8) ) + abs(1-new_alpha-1e-8) * log( matrix(abs(1-new_pai-1e-8), n, p, byrow = TRUE) / abs(1-new_alpha-1e-8) )
      fun3 = new_alpha * ( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log1p(exp(e_mat)) ) + (1-new_alpha)*prob_uncertainty(m, R, matrix(new_zeta, n, p, byrow = TRUE), islog=TRUE)
      y = fun1 + sum(fun2) + sum(fun3)
      y
    }
    
    log_base = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
      new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
      new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
      new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
      model_coefs = model_coefs[-(1:(p*n_factors))]
      new_B = model_coefs[1:p]
      model_coefs = model_coefs[-(1:p)]
      new_beta = matrix(model_coefs[1:(p*q)], p, q)
      
      new_zeta = zeta
      new_pai = pai
      new_alpha = alpha
      
      ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
      e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
      
      y1 = new_alpha * ( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log1p(exp(ll)) ) + (1-new_alpha)*prob_uncertainty(m, R, matrix(new_alpha, n, p, byrow = TRUE), islog=TRUE)
      sum(y1)
    }
    
    log_base2 = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
      new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
      new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
      new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
      model_coefs = model_coefs[-(1:(p*n_factors))]
      new_B = model_coefs[1:p]
      model_coefs = model_coefs[-(1:p)]
      new_beta = matrix(model_coefs[1:(p*q)], p, q)
      
      new_zeta = zeta
      new_pai = pai
      new_alpha = alpha
      
      ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
      e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
      
      y1 = new_alpha * ( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log1p(exp(e_mat)) ) + (1-new_alpha)*prob_uncertainty(m, R, matrix(new_alpha, n, p, byrow = TRUE), islog=TRUE)
      sum(y1)
    }
    
    # update pai
    new_pai = round(apply(new_alpha, 2, mean), 6)
    # update model_coefs: B, Lambda
    fn_model_coefs = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
      new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
      new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
      new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
      model_coefs = model_coefs[-(1:(p*n_factors))]
      new_B = model_coefs[1:p]
      model_coefs = model_coefs[-(1:p)]
      new_beta = matrix(model_coefs[1:(p*q)], p, q)
      
      new_zeta = zeta
      new_pai = pai
      new_alpha = alpha
      
      ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
      e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
      y1 = new_alpha * (m-R) * ll
      y2 = -(m-1)*new_alpha*log1p(exp(e_mat))
      sum(y1) + sum(y2)
    }
    
    gr_model_coefs = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
      new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
      new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
      new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
      model_coefs = model_coefs[-(1:(p*n_factors))]
      new_B = model_coefs[1:p]
      model_coefs = model_coefs[-(1:p)]
      new_beta = matrix(model_coefs[1:(p*q)], p, q)
      
      new_zeta = zeta
      new_pai = pai
      new_alpha = alpha
      
      ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
      e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
      b1 = new_alpha * (m-R) - (m-1)*new_alpha*plogis(e_mat)
      b2 = NULL
      for (w in 1:n_factors) {
        b2 = c(b2, new_alpha*( sweep(m-R, 1, new_mu[, w], "*") - (m-1)*sweep((new_sigma[, w]%*%t(new_Lambda[, w])), 1, new_mu[, w], "+")*plogis(e_mat) ))
      }
      b2 = matrix(b2, n, p*n_factors)
      b3 = NULL
      for (w in 1:q) {
        b3 = c(b3, new_alpha*( sweep(m-R - (m-1)*plogis(e_mat), 1, x1[, w], "*") ))
      }
      b3 = matrix(b3, n, p*q)
      return(c(colSums(b2), colSums(b1), colSums(b3)))
    }
    
    qq = try(optim(c(Lambda, B, beta), fn = fn_model_coefs, gr = gr_model_coefs, method = "BFGS", va_mu = c(new_mu), va_sigma = c(new_sigma), pai = new_pai, alpha = new_alpha, zeta = new_zeta, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
    if ("try-error" %in% class(qq)) {
      new_Lambda = Lambda
      new_B = B
      new_beta = beta
    } else {
      if (iter > 1 && b.cur.logfunc > qq$value) {
        if (trace) {cat("Optimization of model coefs did not improve on iteration step ", iter, "\n")}
        new_Lambda = Lambda
        new_B = B
        new_beta = beta
      } else {
        if (trace) {cat("Model parameters updated", "\n")}
        new_Lambda = matrix(qq$par[1:(p*n_factors)], p, n_factors)
        qq$par = qq$par[-(1:(p*n_factors))]
        new_B = qq$par[1:p]
        qq$par = qq$par[-(1:p)]
        new_beta = matrix(qq$par[1:(p*q)], p, q)
        if (qq$convergence != 0) {
          if (trace) {
            cat("Optimization of model coefs did not converge on iteration step ", iter, "\n")
          }
        }
      }
    }
    
    fn_zeta = function(x, alpha, idx) {
      -sum( (1-alpha[, idx]) * prob_uncertainty(m, R[, idx], x, islog = TRUE) )
    }
    for (j in 1:p) {
      new_zeta[j] = optim(zeta[j], fn_zeta, alpha = new_alpha, idx = j, method = "BFGS", lower = 1e-6)$par
    }
    
    # fn_zeta = function(x, alpha) {
    #   -sum( (1-alpha)*prob_uncertainty(m, R, matrix(x, n, p, byrow = TRUE), islog=TRUE) )
    # }
    # 
    # new_zeta = optim(zeta, fn_zeta, alpha = new_alpha, method = "L-BFGS-B", lower = rep(1e-6, p))$par
    
    # update variational parameter
    delta.alpha.required = 1e-3; p.iter = 1; p.max.iter = 10; delta.alpha = abs(alpha)
    while (!all(delta.alpha < delta.alpha.required) & (p.iter < p.max.iter)) {
      # update alpha
      ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
      e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
      pai_mat = matrix(new_pai, n, p, byrow = TRUE)
      new_alpha = plogis( log(pai_mat / (1-pai_mat)) + log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log1p(exp(e_mat)) - prob_uncertainty(m, R, matrix(new_zeta, n, p, byrow = TRUE), islog=TRUE))
      
      # update sigma
      fn_va_sigma = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
        new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
        new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
        new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
        model_coefs = model_coefs[-(1:(p*n_factors))]
        new_B = model_coefs[1:p]
        model_coefs = model_coefs[-(1:p)]
        new_beta = matrix(model_coefs[1:(p*q)], p, q)
        
        new_zeta = zeta
        new_pai = pai
        new_alpha = alpha
        
        ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
        e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
        y1 = -0.5*( sum(new_sigma) - sum(log(new_sigma)) )
        y2 = -(m-1)*new_alpha*log1p(exp(e_mat))
        y = y1 + sum(y2)
        return(y)
      }
      
      gr_va_sigma = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
        new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
        new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
        new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
        model_coefs = model_coefs[-(1:(p*n_factors))]
        new_B = model_coefs[1:p]
        model_coefs = model_coefs[-(1:p)]
        new_beta = matrix(model_coefs[1:(p*q)], p, q)
        
        new_zeta = zeta
        new_pai = pai
        new_alpha = alpha
        
        ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
        e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
        beta2 = new_Lambda^2
        mu_mat = -(m-1)*new_alpha*plogis(e_mat)
        
        grad.sigma = matrix(NA, n, n_factors)
        if (n_factors == 1) {
          grad.sigma[i, ] = -0.5*(1 - new_sigma[i, ]^(-1)) + apply(mu_mat[i, ]*beta2, 2, sum)
        } else {
          for (i in 1:n) {
            grad.sigma[i, ] = diag(-0.5*(diag(rep(1, n_factors)) - (diag(new_sigma[i, ]))^(-1)) + diag(apply(mu_mat[i, ]*beta2, 2, sum)))
          }
        }
        return(c(grad.sigma))
      }
      
      qq = try(constrOptim(c(sigma), method = "BFGS", f = fn_va_sigma, gr = gr_va_sigma, model_coefs = c(new_Lambda, new_B, new_beta), va_mu = c(new_mu), pai = new_pai, alpha = new_alpha, zeta = new_zeta, ui = diag(1, n*n_factors, n*n_factors), ci = rep(1e-8, n*n_factors), control = list(trace = 0, fnscale = -1, maxit = maxit, reltol = 1e-6)), silent = TRUE)
      if ("try-error" %in% class(qq)) {
        new_sigma = sigma
      } else {
        if (iter > 1 && s.cur.logfunc > qq$value) {
          if (trace) {cat("Optimization of sigma did not improve on iteration step ", iter, "\n")}
          new_sigma = sigma
        } else {
          if (trace) {cat("Variational parameters sigma updated", "\n")}
          new_sigma = matrix(qq$par, n, n_factors)
          if (qq$convergence != 0) {
            if (trace) {
              cat("Optimization of sigma did not converge on iteration step ", iter, "\n")
            }
          }
        }
      }
      
      # update mu
      fn_va_mu = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
        new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
        new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
        new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
        model_coefs = model_coefs[-(1:(p*n_factors))]
        new_B = model_coefs[1:p]
        model_coefs = model_coefs[-(1:p)]
        new_beta = matrix(model_coefs[1:(p*q)], p, q)
        
        new_zeta = zeta
        new_pai = pai
        new_alpha = alpha
        
        ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
        e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
        
        fun1 = -0.5* sum(new_mu^2)
        fun2 = new_alpha * ( (m-R)*new_mu %*% t(new_Lambda) - (m-1)*log1p(exp(e_mat)) )
        y = fun1 + sum(fun2) 
        y
      }
      
      gr_va_mu = function(model_coefs, va_mu, va_sigma, pai, alpha, zeta) {
        new_mu = matrix(va_mu, n, n_factors) # specified to case when n_factors=1
        new_sigma = matrix(va_sigma, n, n_factors) # specified to case when n_factors=1
        new_Lambda = matrix(model_coefs[1:(p*n_factors)], p, n_factors)
        model_coefs = model_coefs[-(1:(p*n_factors))]
        new_B = model_coefs[1:p]
        model_coefs = model_coefs[-(1:p)]
        new_beta = matrix(model_coefs[1:(p*q)], p, q)
        
        new_zeta = zeta
        new_pai = pai
        new_alpha = alpha
        
        ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
        e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
        grad.m = NULL
        sum1 = new_alpha*( (m-R)-(m-1)*plogis(e_mat) )
        for (w in 1:n_factors) {
          grad.m = c(grad.m, rowSums(sweep(sum1, 2, new_Lambda[, w], "*")) - new_mu[, w])
        }
        return(c(grad.m))
      }
      
      qq = try(optim(c(mu), method = "BFGS", fn = fn_va_mu, gr = gr_va_mu, model_coefs = c(new_Lambda, new_B, new_beta), va_sigma = c(new_sigma), pai = new_pai, alpha = new_alpha, zeta = new_zeta, control = list(trace = 0, fnscale = -1, maxit = maxit, reltol = 1e-6)), silent = TRUE)
      if ("try-error" %in% class(qq)) {
        new_mu = mu
      } else {
        if (iter > 1 && m.cur.logfunc > qq$value) {
          if (trace) {cat("Optimization of mu did not improve on iteration step ", iter, "\n")}
          new_mu = mu
        } else {
          if (trace) {cat("Variational parameters mu updated", "\n")}
          new_mu = matrix(qq$par, n, n_factors)
          if (qq$convergence != 0) {
            if (trace) {
              cat("Optimization of mu did not converge on iteration step ", iter, "\n")
            }
          }
        }
      }
      
      q_m = list(value = fn_va_mu(c(new_Lambda, new_B, new_beta), va_mu = c(new_mu), va_sigma = c(new_sigma), pai = new_pai, alpha = new_alpha, zeta = new_zeta))
      m.new.logfunc = q_m$value
      m.cur.logfunc = m.new.logfunc
      
      q_s = list(value = fn_va_sigma(c(new_Lambda, new_B, new_beta), va_mu = c(new_mu), va_sigma = c(new_sigma), pai = new_pai, alpha = new_alpha, zeta = new_zeta))
      s.new.logfunc = q_s$value
      s.cur.logfunc = s.new.logfunc
      
      delta.alpha = abs(new_alpha - alpha)
      sigma = new_sigma
      mu = new_mu
      alpha = new_alpha
      pai = new_pai
      p.iter = p.iter + 1
      
    }
    
    q_b = list(value = fn_model_coefs(c(new_Lambda, new_B, new_beta), va_mu = c(new_mu), va_sigma = c(new_sigma), pai = new_pai, alpha = new_alpha, zeta = new_zeta))
    b.new.logfunc = q_b$value
    b.cur.logfunc = b.new.logfunc
    
    # Take values of VLB to define stopping rule
    qq = list(value = VLB(c(new_Lambda, new_B, new_beta), va_mu = c(new_mu), va_sigma = c(new_sigma), pai = new_pai, alpha = new_alpha, zeta = new_zeta))
    new.VLB = qq$value
    diff = abs(new.VLB - cur.VLB)
    ratio = abs(new.VLB / cur.VLB)
    if (trace) cat("New VLB: ", new.VLB, "cur VLB: ", cur.VLB, "Ratio of VLB: ", ratio, ". Difference in VLB: ", diff, "\n")
    if (trace) cat("\nThe current evaluate results: ", evaluate_Ex4(promax(new.Lambda)$loadings,promax(data$Lambda)$loadings))
    cur.VLB = new.VLB
    
    qq = list(value = log_base(c(new_Lambda, new_B, new_beta), va_mu = c(new_mu), va_sigma = c(new_sigma), pai = new_pai, alpha = new_alpha, zeta = new_zeta))
    cur.log = qq$value
    
    EE = log_base2(c(new_Lambda, new_B, new_beta), va_mu = c(new_mu), va_sigma = c(new_sigma), pai = new_pai, alpha = new_alpha, zeta = new_zeta)
    EN2 = 2*cur.log - 2*EE
    
    
    pai = new_pai
    alpha = new_alpha
    B = new_B
    Lambda = new_Lambda
    mu = new_mu
    sigma = new_sigma
    zeta = new_zeta
    
    iter = iter + 1
  }
  if (iter > 99) {
    print("FACAUB not converging!")
  }
  
  ll = matrix(new_B, n, p, byrow = TRUE) + x1 %*% t(new_beta) + new_mu %*% t(new_Lambda)
  e_mat = ll + 0.5*new_sigma %*% t(new_Lambda^2)
  
  if (is.null(x1)) beta = NULL
  ll = plogis(ll)
  e_mat = plogis(e_mat)
  
  out.list = list()
  out.list$VLB = cur.VLB
  out.list$EE = EE
  out.list$lob = cur.log
  out.list$EN2 = EN2
  out.list$iter = iter - 1
  out.list$lvs$alpha = new_alpha
  out.list$lvs$mu = new_mu
  out.list$lvs$sigma = new_sigma
  out.list$params$pai = new_pai
  out.list$params$B = new_B
  out.list$params$Lambda = new_Lambda
  out.list$params$beta = new_beta
  out.list$params$zeta = new_zeta
  out.list$ll = ll
  out.list$e_mat = e_mat
  
  return(out.list)
}

#Vector and trace correlation coefficients (Ye and Weiss, JASA, 2003)
eval.space <- function(A, B, orthnm = TRUE) 
{
  if(!is.matrix(A)) A <- as.matrix(A)
  if(!is.matrix(B)) B <- as.matrix(B)
  
  if(orthnm)
  { 
    A <- qr.Q(qr(A))
    B <- qr.Q(qr(B)) 
  }
  
  mat <- t(B) %*% A %*% t(A) %*% B
  d <- eigen(mat)$values
  d <- (d + abs(d))/2
  q <- sqrt(prod(d))
  r <- sqrt(mean(d))
  
  c(q, r, acos(q))
}


generate_data_FACUB = function(n, p, q, d, m, b0, pai_param, type1="uniform", type2 = "identity") {
  a = rep(0, n)  # 个体固定效应a_i, i=1,\dots,n
  b = rep(b0, p)  # 变量的固定效应b_j, j=1,\dots,p
  b = rnorm(p)
  Lambda = matrix(NA, p, d)  # 因子载荷矩阵Lambda pxd维矩阵，从均匀分布或标准正态中生成
  if (type1 == "uniform") {
    for (i in 1:p) {
      Lambda[i, ] = runif(d, -2, 2)
    }
  } else if (type1 == "gaussian") {
    for (i in 1:p) {
      Lambda[i, ] = rnorm(d, 0, 1)
    }
  }
  
  f = matrix(NA, n, d)  # 因子f nxd维矩阵，从标准正态中生成
  for (i in 1:n) {
    f[i, ] = rnorm(d, 0, 1)
  }
  if (q == 0) {   # 没有协变量
    x = NULL
    beta = NULL
    csi = plogis(matrix(a, n, p) + matrix(b, n, p, byrow = TRUE) + f %*% t(Lambda))
  } else {
    beta = matrix(rnorm(p*q), p, q)
    if (type2 == "identity") {
      x = matrix(rnorm(n*q), n, q)  
    } else if (type2 == "exchange") {
      rho = .3
      Sigma = matrix(rho, q, q)
      diag(Sigma) = rep(1, q)
      x = mvrnorm(n, mu = rep(0, q), Sigma = Sigma)
    }
    csi = plogis(matrix(a, n, p) + matrix(b, n, p, byrow = TRUE) + f %*% t(Lambda) + x %*% t(beta))
  }
  
  pai = runif(p, pai_param[1], pai_param[2])   # CUB的参数pai
  z = R = matrix(NA, n, p)
  for (j in 1:p) {
    z[, j] = sample(0:1, n, replace = TRUE, prob = c(1-pai[j], pai[j]))
    z_idx = which(z[, j] == 1)
    R[z_idx, j] = rbinom(length(z_idx), m-1, prob = 1-csi[z_idx, j]) + 1
    z0_idx = setdiff(1:n, z_idx)
    R[z0_idx, j] = sample(1:m, length(z0_idx), replace = TRUE, prob = rep(1/m, m))
  }
  pai_mat = matrix(pai, n, p, byrow = TRUE)
  prob = pai_mat*prob_feeling(m, R, csi) + (1-pai_mat)/m
  # nk = numeric(n)
  # for (i in 1:n) {
  #   z[i, ] = rbinom(p, size = 1, prob = pai)
  #   nk[i] = length(which(z[i, ] == 0))
  #   for (j in 1:p) {
  #     R[i, j] = rbinom(1, size = m-1, prob = 1 - csi[i, j]) + 1
  #   }
  #   R[i, z[i, ] == 0] = sample(1:m, nk[i], replace = TRUE, prob = rep(1/m, m))
  # }
  
  return(list(R = R, x = x, beta = beta, pai = pai, Lambda = Lambda, f = f, A = a, B = b, csi = csi, prob = prob))
}


generate_data_FACUB2 = function(n, p, q, d, m, b0, pai_param1, pai_param2, type1="uniform", type2 = "identity") {
  a = rep(0, n)  # 个体固定效应a_i, i=1,\dots,n
  b = rep(b0, p)  # 变量的固定效应b_j, j=1,\dots,p
  b = rnorm(p)
  Lambda = matrix(NA, p, d)  # 因子载荷矩阵Lambda pxd维矩阵，从均匀分布或标准正态中生成
  if (type1 == "uniform") {
    for (i in 1:p) {
      Lambda[i, ] = runif(d, -2, 2)
    }
  } else if (type1 == "gaussian") {
    for (i in 1:p) {
      Lambda[i, ] = rnorm(d, 0, 1)
    }
  }
  
  f = matrix(NA, n, d)  # 因子f nxd维矩阵，从标准正态中生成
  for (i in 1:n) {
    f[i, ] = rnorm(d, 0, 1)
  }
  if (q == 0) {   # 没有协变量
    x = NULL
    beta = NULL
    csi = plogis(matrix(a, n, p) + matrix(b, n, p, byrow = TRUE) + f %*% t(Lambda))
  } else {
    beta = matrix(rnorm(p*q), p, q)
    if (type2 == "identity") {
      x = matrix(rnorm(n*q), n, q)  
    } else if (type2 == "exchange") {
      rho = .3
      Sigma = matrix(rho, q, q)
      diag(Sigma) = rep(1, q)
      x = mvrnorm(n, mu = rep(0, q), Sigma = Sigma)
    }
    csi = plogis(matrix(a, n, p) + matrix(b, n, p, byrow = TRUE) + f %*% t(Lambda) + x %*% t(beta))
  }
  
  # z2 = sample(0:1, p, replace=TRUE, prob=c(1/3,2/3))
  # pai = numeric(p)
  # pai[z2 == 0] = runif(length(which(z2 == 0)), pai_param1[1], pai_param1[2])
  # pai[z2 == 1] = runif(length(which(z2 == 1)), pai_param2[1], pai_param2[2])
  pai1 = runif(floor(p/3), pai_param1[1], pai_param1[2])   # CUB的参数pai
  pai2 = runif(p-length(pai1), pai_param2[1], pai_param2[2])
  pai = c(pai1, pai2)
  z = R = matrix(NA, n, p)
  for (j in 1:p) {
    z[, j] = sample(0:1, n, replace = TRUE, prob = c(1-pai[j], pai[j]))
    z_idx = which(z[, j] == 1)
    R[z_idx, j] = rbinom(length(z_idx), m-1, prob = 1-csi[z_idx, j]) + 1
    z0_idx = setdiff(1:n, z_idx)
    R[z0_idx, j] = sample(1:m, length(z0_idx), replace = TRUE, prob = rep(1/m, m))
  }
  pai_mat = matrix(pai, n, p, byrow = TRUE)
  prob = pai_mat*prob_feeling(m, R, csi) + (1-pai_mat)/m
  # nk = numeric(n)
  # for (i in 1:n) {
  #   z[i, ] = rbinom(p, size = 1, prob = pai)
  #   nk[i] = length(which(z[i, ] == 0))
  #   for (j in 1:p) {
  #     R[i, j] = rbinom(1, size = m-1, prob = 1 - csi[i, j]) + 1
  #   }
  #   R[i, z[i, ] == 0] = sample(1:m, nk[i], replace = TRUE, prob = rep(1/m, m))
  # }
  
  return(list(R = R, x = x, beta = beta, pai = pai, Lambda = Lambda, f = f, A = a, B = b, csi = csi, prob = prob))
}

generate_initials = function(R, x, n.factors, type="1", B0=1, beta0=0, sigma0=0.01, pai0=0.5) {
  n = nrow(R)
  p = ncol(R)
  q = ifelse(is.null(x), 1, ncol(x))
  
  if (type == "1") {
    R_sc = scale(log2(1+R), center = TRUE, scale = TRUE)
    re = svd(R_sc, n.factors, n.factors)
    mu = re$u %*% diag(re$d[1:n.factors], nrow = n.factors)
    Lambda = re$v
  } else if (type == "2") {
    poly_res = polychoric(R)
    R_poly = poly_res$rho
    re = svd(R_poly, n.factors, n.factors)
    mu = re$u %*% diag(re$d[1:n.factors], nrow = n.factors)
    Lambda = re$v
  } else if (type == "3") {
    R_norm = apply(R, 2, function(x) qnorm(rank(x) / (length(x) + 1)))
    R_scaled = scale(R_norm)
    re = svd(R_scaled, n.factors, n.factors)
    mu = re$u * sqrt(n-1)
    Lambda = re$v %*% diag(re$d[1:n.factors], nrow = n.factors) / sqrt(n-1)
  }
  
  B = rep(B0, p)
  beta = matrix(beta0, p, q)
  sigma = matrix(sigma0, n, n.factors)
  pai = rep(pai0, p)
  alpha = matrix(pai, n, p, byrow = TRUE)
  list(mu = mu, sigma = sigma, pai = pai, alpha = alpha, B = B, beta = beta, Lambda = Lambda)
}

FAVA_new = function(R, V=NULL, m, n.factors, initials, maxit=200, trace=FALSE) {
  #newR = check_R(R, m)
  #R = newR$R
  #re_order = newR$reorder
  n = nrow(R)
  p = ncol(R)
  # if covariates included
  if (is.null(V)) {
    q = 1
    X = matrix(0, n, q)
  } else {
    X = V
    q = ncol(X)
  }
  out.list = list()
  
  mu = new.mu = initials$mu
  sigma = new.sigma = initials$sigma
  pai = new.pai = initials$pai
  alpha = new.alpha = initials$alpha
  B = new.B = initials$B
  beta = new.beta = initials$beta
  Lambda = new.Lambda = initials$Lambda
  
  # VA iteration about variational lower bound
  cur.VLB = -1e6; iter = 1; ratio = 10; diff = 1e5; eps = 1e-4; max.iter = 100;
  m.cur.logfunc = -1e6; b.cur.logfunc = -1e6; d.cur.logfunc = -1e6; s.cur.logfunc = -1e6;
  
  while ( (diff > eps*(abs(cur.VLB)+eps)) && iter <= max.iter ) {
    if (trace) {
      cat("=======================\nIteration: ", iter, "\n")
      cat(" [Outer ", iter, "] j=1 Lambda START: ", new.Lambda[1, ], "\n")
      cat(" [Outer ", iter, "]  pai START: ", new.pai, "\n")
      cat(" [Outer ", iter, "]  B START: ", new.B, "\n")
    }
    
    ## VLB b   
    VLB = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      
      fun1 = -0.5*( sum(new.sigma) + sum(new.mu^2) - sum(log(new.sigma)) )
      fun2 = (new.alpha+1e-8) * log( matrix(new.pai+1e-8, n, p, byrow=TRUE) / (new.alpha+1e-8) ) + abs(1-new.alpha-1e-8) * log( matrix(abs(1-new.pai-1e-8), n, p, byrow=TRUE) / abs(1-new.alpha-1e-8) ) 
      fun3 = new.alpha*( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(e.mat)) ) - (1-new.alpha)*log(m)
      
      y = fun1 + sum(fun2) + sum(fun3)
      return(y)
    }
    
    # log-base
    log_base = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      y1 = new.alpha*( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(ll)) ) - (1-new.alpha)*log(m)
      return(sum(y1))
    }
    
    # log-base
    log_base2 = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      y1 = new.alpha*( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(e.mat)) ) - (1-new.alpha)*log(m)
      return(sum(y1))
    }
    
    # update pai
    new.pai = round(apply(new.alpha, 2, mean), 6)
    
    # update model.coefs: B, beta, Lambda
    model.coefs_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      y1 = new.alpha*(m-R)*ll
      y2 = -new.alpha*(m-1)*log(1+exp(e.mat))
      y = sum(y1) + sum(y2)
      return(y)
    }
    
    model.coefs_gr_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      b1 = new.alpha*(m-R) - (m-1)*new.alpha*plogis(e.mat)
      b2 = NULL
      for (w in 1:n.factors) {
        b2 = c(b2, new.alpha*( sweep(m-R, 1, new.mu[, w], "*") - (m-1)*sweep((new.sigma[, w]%*%t(new.Lambda[, w])), 1, new.mu[, w], "+")*exp(e.mat)/(1+exp(e.mat)) ))
      }
      b2 = matrix(b2, n, p*n.factors)
      b3 = NULL
      for (w in 1:q) {
        b3 = c(b3, new.alpha*( sweep(m-R - (m-1)*plogis(e.mat), 1, X[, w], "*" )))
      }
      b3 = matrix(b3, n, p*q)
      return(c(colSums(b2), colSums(b1), colSums(b3)))
    }
    
    qq = try(optim(c(Lambda, B, beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha, method = "BFGS", fn = model.coefs_f, gr = model.coefs_gr_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
    
    if ("try-error" %in% class(qq)) {
      new.Lambda = Lambda; new.B = B; new.beta = beta
    } else {
      if (iter > 1 && b.cur.logfunc > qq$value) {
        if (trace) {cat("Optimization of model coefs dit not improve on iteration step ", iter, "\n")}
        new.Lambda = Lambda; new.B = B; new.beta = beta
      } else {
        if (trace) {cat("Model parameters updated", "\n")}
        new.Lambda = matrix(qq$par[1:(p*n.factors)], p, n.factors); qq$par = qq$par[-(1:(p*n.factors))]
        new.B = qq$par[1:p]; qq$par = qq$par[-(1:p)]
        new.beta = matrix(qq$par[1:(p*q)], p, q)
        if (qq$convergence != 0) {
          if (trace) {
            cat("Optimization of model coefs did not converge on iteration step ", iter, "\n")
          }
        }
      }
    }
    
    delta.alpha.required = 1e-3; p.iter = 1; p.max.iter = 10; delta.alpha = abs(alpha)
    while (!all(delta.alpha < delta.alpha.required) & (p.iter < p.max.iter)) {
      # update alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      pai.mat = matrix(new.pai, n, p, byrow = TRUE)
      new.alpha = plogis( log(pai.mat / (1-pai.mat)) + log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(e.mat)) + log(m) )
      if (trace) cat(" [Inner ", p.iter, "] alpha START: ", new.alpha[1, ], "\n")
      
      # update sigma
      va.sigma_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
        new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p*q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
        
        fun1 = -0.5*( sum(new.sigma) - sum(log(new.sigma)) )
        y2 = -(m-1)*new.alpha*log(1+exp(e.mat))
        y = fun1 + sum(y2)
        return(y)
      }
      
      va.sigma_gr_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
        new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p*q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
        beta2 = new.Lambda^2
        mu.mat = -(m-1)*new.alpha*exp(e.mat)/(1+exp(e.mat))
        
        grad.sigma = matrix(NA, n, n.factors)
        if (n.factors == 1) {
          for (i in 1:n) {
            grad.sigma[i, ] = -0.5*(1 - new.sigma[i, ]^(-1)) + apply(mu.mat[i, ]*beta2, 2, sum)
          }
        } else {
          for (i in 1:n) {
            grad.sigma[i, ] = diag(-0.5*(diag(rep(1, n.factors)) - diag(1/new.sigma[i, ], n.factors)) + diag(apply(mu.mat[i, ]*beta2, 2, sum)))
          }
        }
        return(c(grad.sigma))
      }
      if (trace) cat(" [Inner ", p.iter, "] i=0 SIGMA START: ", new.sigma[1, ], "\n")
      qq = try(constrOptim(c(sigma), method = "BFGS", f = va.sigma_f, gr = va.sigma_gr_f, model.coefs = c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), pai = new.pai, alpha = new.alpha, ui = diag(1, n*n.factors, n*n.factors), ci = rep(1e-8, n*n.factors), control = list(trace = 0, fnscale = -1, maxit = maxit, reltol = 1e-6)), silent = TRUE)
      
      if ("try-error" %in% class(qq)) {
        new.sigma = sigma
      } else {
        if (iter > 1 && s.cur.logfunc > qq$value) {
          if (trace) {cat("Optimization of sigma did not improve on iteration step ", iter, "\n")}
          new.sigma = sigma
        } else {
          if (trace) {cat("Variational parameters sigma updated", "\n")}
          new.sigma = matrix(qq$par, n, n.factors)
          if (qq$convergence != 0) {
            if (trace) {
              cat("Optimization of sigma did not converge on iteration step ", iter, "\n")
            }
          }
        }
      }
      
      # update mu
      mu_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
        new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p*q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
        fun1 = -0.5*( sum(new.mu^2) )
        fun2 = new.alpha * ( (m-R)*new.mu %*% t(new.Lambda) - (m-1)*log(1+exp(e.mat)) )
        y = fun1 + sum(fun2)
        return(y)
      }
      
      mu_grad_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
        new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p*q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
        grad.m = NULL
        sum1 = new.alpha*( (m-R)-(m-1)*exp(e.mat)/(1+exp(e.mat)) )
        for (w in 1:n.factors) {
          grad.m = c(grad.m, rowSums(sweep(sum1, 2, new.Lambda[, w], "*")) - new.mu[, w])
        }
        return(c(grad.m))
      }
      if (trace) cat(" [Inner ", p.iter, "] i=0 MU START: ", new.mu[1, ], "\n")
      #qq = try(optim(c(mu), method = "BFGS", fn = mu_f, gr = mu_grad_f, model.coefs = c(new.Lambda, new.B, new.beta), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha, control = list(trace = 0, fnscale = -1, maxit = maxit, reltol = 1e-6)), silent = TRUE)
      qq = try(optim(c(mu), method = "BFGS", fn = mu_f, gr = mu_grad_f, model.coefs = c(new.Lambda, new.B, new.beta), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha, control = list(trace = 0, fnscale = -1, maxit = maxit, reltol = 1e-6)), silent = TRUE)
      
      if ("try-error" %in% class(qq)) {
        new.mu = mu
      } else {
        if (iter > 1 && m.cur.logfunc > qq$value) {
          if (trace) {cat("Optimization of mu did not improve on iteration step ", iter, "\n")}
          new.mu = mu
        } else {
          if (trace) {cat("Variational parameters mu updated", "\n")}
          new.mu = matrix(qq$par, n, n.factors)
          if (qq$convergence != 0) {
            if (trace) {
              cat("Optimization of mu did not converge on iteration step ", iter, "\n")
            }
          }
        }
      }
      
      q_m = list(value = mu_f(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
      m.new.logfunc = q_m$value
      m.cur.logfunc = m.new.logfunc
      
      q_s = list(value = va.sigma_f(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
      s.new.logfunc = q_s$value
      s.cur.logfunc = s.new.logfunc
      
      delta.alpha = abs(new.alpha - alpha)
      sigma = new.sigma
      mu = new.mu
      alpha = new.alpha
      pai = new.pai
      p.iter = p.iter + 1
    }
    
    q_b = list(value = model.coefs_f(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
    b.new.logfunc = q_b$value
    b.cur.logfunc = b.new.logfunc
    
    # Take values of VLB to define stopping rule
    qq = list(value = VLB(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
    new.VLB = qq$value
    diff = abs(new.VLB - cur.VLB)
    ratio = abs(new.VLB / cur.VLB)
    if (trace) cat("New VLB: ", new.VLB, "cur VLB: ", cur.VLB, "Ratio of VLB: ", ratio, ". Difference in VLB: ", diff, "\n")
    #if (trace) cat("\nThe current evaluate results: ", evaluate_Ex4(promax(new.Lambda)$loadings, promax(data$Lambda)$loadings))
    cur.VLB = new.VLB
    
    qq = list(value = log_base(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
    cur.log = qq$value
    
    EE = log_base2(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha)
    EN2 = 2*cur.log - 2*EE
    
    
    pai = new.pai
    alpha = new.alpha
    B = new.B
    Lambda = new.Lambda
    beta = new.beta
    mu = new.mu
    sigma = new.sigma
    
    iter = iter + 1
  }
  
  if (iter > 99) {
    print("FACUB not converging!")
  }
  
  
  ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
  e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
  
  
  ll = plogis(ll)
  e.mat = plogis(e.mat)
  #if (re_order) {
  #  B = -B
  #  beta = -beta
  #  e.mat = 1-e.mat
  #  ll = 1-ll
  #}
  
  if (is.null(V)) beta = NULL
  
  out.list$VLB = cur.VLB
  out.list$EE = EE
  out.list$lob = cur.log
  out.list$EN2 = EN2
  out.list$iter = iter - 1
  out.list$lvs$alpha = alpha
  out.list$lvs$mu = mu
  out.list$lvs$sigma = sigma
  out.list$params$pai = pai
  out.list$params$B = B
  out.list$params$beta = beta
  out.list$params$Lambda = Lambda
  out.list$ll = ll
  out.list$e.mat = e.mat
  
  return(out.list)
}

FAVA_new2 = function(R, V=NULL, m, n.factors, initials, maxit=200, trace=FALSE) {
  #newR = check_R(R, m)
  #R = newR$R
  #re_order = newR$reorder
  n = nrow(R)
  p = ncol(R)
  # if covariates included
  if (is.null(V)) {
    q = 1
    X = matrix(0, n, q)
  } else {
    X = V
    q = ncol(X)
  }
  out.list = list()
  
  mu = new.mu = initials$mu
  sigma = new.sigma = initials$sigma
  pai = new.pai = initials$pai
  alpha = new.alpha = initials$alpha
  B = new.B = initials$B
  beta = new.beta = initials$beta
  Lambda = new.Lambda = initials$Lambda
  
  # VA iteration about variational lower bound
  cur.VLB = -1e6; iter = 1; ratio = 10; diff = 1e5; eps = 1e-4; max.iter = 100;
  m.cur.logfunc = -1e6; b.cur.logfunc = -1e6; d.cur.logfunc = -1e6; s.cur.logfunc = -1e6;
  
  while ( (diff > eps*(abs(cur.VLB)+eps)) && iter <= max.iter ) {
    if (trace) {
      cat("=======================\nIteration: ", iter, "\n")
      cat(" [Outer ", iter, "] j=1 Lambda START: ", new.Lambda[1, ], "\n")
      cat(" [Outer ", iter, "]  pai START: ", new.pai, "\n")
      cat(" [Outer ", iter, "]  B START: ", new.B, "\n")
    }
    ## VLB b   
    VLB = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      
      fun1 = -0.5*( sum(new.sigma) + sum(new.mu^2) - sum(log(new.sigma)) )
      fun2 = (new.alpha+1e-8) * log( matrix(new.pai+1e-8, n, p, byrow=TRUE) / (new.alpha+1e-8) ) + abs(1-new.alpha-1e-8) * log( matrix(abs(1-new.pai-1e-8), n, p, byrow=TRUE) / abs(1-new.alpha-1e-8) ) 
      fun3 = new.alpha*( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(e.mat)) ) - (1-new.alpha)*log(m)
      
      y = fun1 + sum(fun2) + sum(fun3)
      return(y)
    }
    
    # log-base
    log_base = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      y1 = new.alpha*( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(ll)) ) - (1-new.alpha)*log(m)
      return(sum(y1))
    }
    
    # log-base
    log_base2 = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      y1 = new.alpha*( log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(e.mat)) ) - (1-new.alpha)*log(m)
      return(sum(y1))
    }
    
    # update pai
    new.pai = round(apply(new.alpha, 2, mean), 6)
    
    # update model.coefs: B, beta, Lambda
    model.coefs_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      y1 = new.alpha*(m-R)*ll
      y2 = -new.alpha*(m-1)*log(1+exp(e.mat))
      y = sum(y1) + sum(y2)
      return(y)
    }
    
    model.coefs_gr_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
      new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p*q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      b1 = new.alpha*(m-R) - (m-1)*new.alpha*exp(e.mat)/(1+exp(e.mat))
      b2 = NULL
      for (w in 1:n.factors) {
        b2 = c(b2, new.alpha*( sweep(m-R, 1, new.mu[, w], "*") - (m-1)*sweep((new.sigma[, w]%*%t(new.Lambda[, w])), 1, new.mu[, w], "+")*exp(e.mat)/(1+exp(e.mat)) ))
      }
      b2 = matrix(b2, n, p*n.factors)
      b3 = NULL
      for (w in 1:q) {
        b3 = c(b3, new.alpha*( sweep(m-R - (m-1)*exp(e.mat)/(1+exp(e.mat)), 1, X[, w], "*" )))
      }
      b3 = matrix(b3, n, p*q)
      return(c(colSums(b2), colSums(b1), colSums(b3)))
    }
    
    qq = try(optim(c(Lambda, B, beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha, method = "BFGS", fn = model.coefs_f, gr = model.coefs_gr_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
    
    if ("try-error" %in% class(qq)) {
      new.Lambda = Lambda; new.B = B; new.beta = beta
    } else {
      if (iter > 1 && b.cur.logfunc > qq$value) {
        if (trace) {cat("Optimization of model coefs dit not improve on iteration step ", iter, "\n")}
        new.Lambda = Lambda; new.B = B; new.beta = beta
      } else {
        #cat("Model parameters updated\n")
        if (trace) {cat("Model parameters updated", "\n")}
        new.Lambda = matrix(qq$par[1:(p*n.factors)], p, n.factors); qq$par = qq$par[-(1:(p*n.factors))]
        new.B = qq$par[1:p]; qq$par = qq$par[-(1:p)]
        new.beta = matrix(qq$par[1:(p*q)], p, q)
        if (qq$convergence != 0) {
          if (trace) {
            cat("Optimization of model coefs did not converge on iteration step ", iter, "\n")
          }
        }
      }
    }
    
    delta.alpha.required = 1e-3; p.iter = 1; p.max.iter = 10; delta.alpha = abs(alpha)
    while (!all(delta.alpha < delta.alpha.required) & (p.iter < p.max.iter)) {
      # update alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
      e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
      pai.mat = matrix(new.pai, n, p, byrow = TRUE)
      new.alpha = plogis( log(pai.mat / (1-pai.mat)) + log(choose(m-1, R-1)) + (m-R)*ll - (m-1)*log(1+exp(e.mat)) + log(m) )
      if (trace) cat(" [Inner ", p.iter, "] alpha START: ", new.alpha[1, ], "\n")
      
      # update sigma
      va.sigma_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
        new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p*q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
        
        fun1 = -0.5*( sum(new.sigma) - sum(log(new.sigma)) )
        y2 = -(m-1)*new.alpha*log(1+exp(e.mat))
        y = fun1 + sum(y2)
        return(y)
      }
      
      va.sigma_gr_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
        new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p*q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
        beta2 = new.Lambda^2
        mu.mat = -(m-1)*new.alpha*plogis(e.mat)
        
        grad.sigma = matrix(NA, n, n.factors)
        if (n.factors == 1) {
          for (i in 1:n) {
            grad.sigma[i, ] = -0.5*(1 - new.sigma[i, ]^(-1)) + apply(mu.mat[i, ]*beta2, 2, sum)
          }
        } else {
          for (i in 1:n) {
            grad.sigma[i, ] = diag(-0.5*(diag(rep(1, n.factors)) - diag(1/new.sigma[i, ], n.factors)) + diag(apply(mu.mat[i, ]*beta2, 2, sum)))
          }
        }
        return(c(grad.sigma))
      }
      
      va.sigma_fi = function(model.coefs, va.mu, va.sigma, 
                             pai, alpha, ii) {
        new.mu = va.mu # n.factors dimension
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = va.sigma # n.factors dimension
        new.pai = pai
        new.alpha = alpha # p dimension
        ll = new.B + X[ii, ] %*% t(new.beta) + new.mu %*% t(new.Lambda) # p dimension
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        y1 = -0.5 * (sum(new.sigma) - sum(log(new.sigma)))
        y2 = -(m - 1) * new.alpha * log(1 + exp(e.mat))
        y = y1 + sum(y2)
        return(y)
      }
      
      va.sigma_gr_fi = function(model.coefs, va.mu, va.sigma, 
                                pai, alpha, ii) {
        new.mu = va.mu # n.factors dimension
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = va.sigma # n.factors dimension
        new.pai = pai
        new.alpha = alpha # p dimension
        ll = c(new.B + X[ii, ] %*% t(new.beta) + new.mu %*% t(new.Lambda)) # p dimension
        e.mat = c(ll + 0.5 * (new.sigma) %*% t(new.Lambda^2))
        beta2 = new.Lambda^2
        mu.vec = -(m - 1) * new.alpha * plogis(e.mat)
        grad.sigma = diag(-0.5 * (diag(rep(1, n.factors), n.factors) - diag(new.sigma^(-1), n.factors)) + 
                            diag(colSums(mu.vec * beta2), n.factors))
        
        return(c(grad.sigma))
      }
      for (i in 1:n) {
        if (i == 1) {
          if (trace) cat(" [Inner ", p.iter, "] i=0 SIGMA START: ", new.sigma[i, ], "\n")
        }
        
        qq = try(constrOptim(sigma[i, ], method = "BFGS", f = va.sigma_fi,
                             gr = va.sigma_gr_fi,
                             model.coefs = c(new.Lambda,new.B, new.beta),
                             va.mu = new.mu[i, ],
                             pai = new.pai,
                             alpha = new.alpha[i, ],
                             ii = i,
                             ui = diag(1, n.factors, n.factors),
                             ci = rep(1e-08, n.factors),
                             control = list(trace = 0,
                                            fnscale = -1,
                                            maxit = maxit,
                                            reltol = 1e-06)), silent = TRUE)
        if ("try-error" %in% class(qq)) {
          new.sigma[i, ] = sigma[i, ]
        } else {
          if (iter > 1 && s.cur.logfunc > qq$value) {
            #if (trace) { cat("Optimization of sigma did not improve on iteration step ", iter, "\n")}
            new.sigma[i, ] = sigma[i, ]
          } else {
            #if (trace) {cat("Variational parameters sigma updated","\n")}
            if (trace) if (i == 1) cat("Variational parameters sigma updated\n")
            new.sigma[i, ] = qq$par
            if (qq$convergence != 0) {
              if (trace) if (i == 1) cat("Optimization of sigma did not converge!\n")
              #if (trace) {cat("Optimization of sigma did not converge on iteration step ",iter, "\n")}
            }
          }
        }
      }
      # update mu
      mu_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
        new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p*q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
        fun1 = -0.5*( sum(new.mu^2) )
        fun2 = new.alpha * ( (m-R)*new.mu %*% t(new.Lambda) - (m-1)*log(1+exp(e.mat)) )
        y = fun1 + sum(fun2)
        return(y)
      }
      
      mu_grad_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p*n.factors)], p, n.factors); model.coefs = model.coefs[-(1:(p*n.factors))]
        new.B = model.coefs[1:p]; model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p*q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
        grad.m = NULL
        sum1 = new.alpha*( (m-R)-(m-1)*exp(e.mat)/(1+exp(e.mat)) )
        for (w in 1:n.factors) {
          grad.m = c(grad.m, rowSums(sweep(sum1, 2, new.Lambda[, w], "*")) - new.mu[, w])
        }
        return(c(grad.m))
      }
      
      mu_fi = function(model.coefs, va.mu, va.sigma, pai, alpha, ii) {
        new.mu = va.mu # n.factors dimension
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = va.sigma # n.factors dimension
        new.pai = pai
        new.alpha = alpha # p dimension
        ll = new.B + X[ii, ] %*% t(new.beta) + new.mu %*% t(new.Lambda) # p dimension
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        y1 = -0.5*sum(new.mu^2)
        y2 = new.alpha * (m - R[ii, ]) * new.mu %*% t(new.Lambda)
        y3 = -(m - 1) * new.alpha * log(1 + exp(e.mat))
        y = y1 + sum(y2) + sum(y3)
        return(y)
      }
      mu_grad_fi = function(model.coefs, va.mu, va.sigma, pai, alpha, ii) {
        new.mu = va.mu # n.factors dimension
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p,q)
        new.sigma = va.sigma # n.factors dimension
        new.pai = pai
        new.alpha = alpha # p dimension
        ll = c(new.B + X[ii, ] %*% t(new.beta) + new.mu %*% t(new.Lambda)) # p dimension
        e.mat = c(ll + 0.5 * (new.sigma) %*% t(new.Lambda^2))
        sum1 = new.alpha * ((m-R[ii, ]) - (m - 1) * plogis(e.mat))
        return(-new.mu + as.vector(crossprod(sum1, new.Lambda)))
      }
      for (i in 1:n) {
        if (i == 1) {
          if (trace) cat(" [Inner ", p.iter, "] i=0 MU START: ", new.mu[i, ], "\n")
        }
        qq = try(optim(mu[i, ], method = "BFGS", fn = mu_fi,
                       gr = mu_grad_fi,
                       model.coefs = c(new.Lambda, new.B, new.beta),
                       va.sigma = new.sigma[i, ],
                       pai = new.pai,
                       alpha = new.alpha[i, ],
                       ii = i,
                       control = list(trace = 0,
                                      fnscale = -1,
                                      maxit = maxit,
                                      reltol = 1e-06)), silent = TRUE)
        if ("try-error" %in% class(qq)) {
          new.mu[i, ] = mu[i, ]
          if (trace) if (i == 1) print("Optimization of mu did not improve on iteration step")
        } else {
          if (iter > 1 && m.cur.logfunc > qq$value) {
            if (trace) if (i == 1) print("Optimization of mu did not improve on iteration step")
            #if (trace) {cat("Optimization of mu did not improve on iteration step ",iter, "\n")}
            new.mu[i, ] = mu[i, ]
          } else {
            if (trace) if (i == 1) cat("Variational parameters mu updated\n")
            new.mu[i, ] = qq$par
            if (qq$convergence != 0) {
              if (trace) if (i == 1) print("Optimization of mu did not improve on iteration step")
            }
          }
        }
      }
      
      
      
      q_m = list(value = mu_f(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
      m.new.logfunc = q_m$value
      m.cur.logfunc = m.new.logfunc
      
      q_s = list(value = va.sigma_f(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
      s.new.logfunc = q_s$value
      s.cur.logfunc = s.new.logfunc
      
      delta.alpha = abs(new.alpha - alpha)
      sigma = new.sigma
      mu = new.mu
      alpha = new.alpha
      pai = new.pai
      p.iter = p.iter + 1
    }
    
    q_b = list(value = model.coefs_f(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
    b.new.logfunc = q_b$value
    b.cur.logfunc = b.new.logfunc
    
    # Take values of VLB to define stopping rule
    qq = list(value = VLB(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
    new.VLB = qq$value
    diff = abs(new.VLB - cur.VLB)
    ratio = abs(new.VLB / cur.VLB)
    if (trace) cat("New VLB: ", new.VLB, "cur VLB: ", cur.VLB, "Ratio of VLB: ", ratio, ". Difference in VLB: ", diff, "\n")
    #if (trace) cat("\nThe current evaluate results: ", evaluate_Ex4(promax(new.Lambda)$loadings, promax(data$Lambda)$loadings))
    cur.VLB = new.VLB
    
    qq = list(value = log_base(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha))
    cur.log = qq$value
    
    EE = log_base2(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha)
    EN2 = 2*cur.log - 2*EE
    
    
    pai = new.pai
    alpha = new.alpha
    B = new.B
    Lambda = new.Lambda
    beta = new.beta
    mu = new.mu
    sigma = new.sigma
    
    iter = iter + 1
  }
  
  if (iter > 99) {
    print("FACUB not converging!")
  }
  
  
  ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
  e.mat = ll + .5* (new.sigma) %*% t(new.Lambda^2)
  
  
  ll = plogis(ll)
  e.mat = plogis(e.mat)
  #if (re_order) {
  #  B = -B
  #  beta = -beta
  #  e.mat = 1-e.mat
  #  ll = 1-ll
  #}
  
  if (is.null(V)) beta = NULL
  
  out.list$VLB = cur.VLB
  out.list$EE = EE
  out.list$lob = cur.log
  out.list$EN2 = EN2
  out.list$iter = iter - 1
  out.list$lvs$alpha = alpha
  out.list$lvs$mu = mu
  out.list$lvs$sigma = sigma
  out.list$params$pai = pai
  out.list$params$B = B
  out.list$params$beta = beta
  out.list$params$Lambda = Lambda
  out.list$ll = ll
  out.list$e.mat = e.mat
  
  return(out.list)
}


error = function(a, b) {
  proj_a = a %*% MASS::ginv(t(a)%*%a) %*% t(a)
  proj_b = b %*% MASS::ginv(t(b)%*%b) %*% t(b)
  error1 = vegan::procrustes(a, b, symmetric = TRUE)$ss   # procrustes error
  error2 = Matrix::norm(proj_a - proj_b, "2")      # projection error
  
  return(c(error1, error2))
}

evaluate_Ex4 = function(est_Lambda, true_Lambda) {
  c(error(est_Lambda, true_Lambda), eval.space(est_Lambda, true_Lambda)[1:2])
}


evaluate_FACUB = function(true, estimated) {
  MSE_pai = MSE(estimated$params$pai, true$pai)
  MSE_B = MSE(estimated$params$B, true$B)
  if (is.null(true$beta)) {
    estimated$params$beta = NULL
    MSE_beta = 0
  } else {
    MSE_beta = MSE(estimated$params$beta, true$beta)
  }
  
  MSE_emat = MSE(estimated$e.mat, true$csi)
  MSE_ll = MSE(estimated$ll, true$csi)
  if (is.null(true$beta)) estimated$params$beta = NULL
  mat1 = cbind(true$B, true$beta, true$Lambda)
  mat2 = cbind(estimated$params$B, estimated$params$beta, estimated$params$Lambda)
  c(evaluate_Ex4(mat1, mat2), MSE_pai, MSE_B, MSE_beta, MSE_emat, MSE_ll)
}

evaluate_FACUB2 = function(true, estimated) {
  MSE_pai = MSE(estimated$params$pai, true$pai)
  MSE_B = MSE(estimated$params$B, true$B)
  MSE_emat = MSE(estimated$e.mat, true$csi)
  MSE_ll = MSE(estimated$ll, true$csi)
  if (is.null(true$beta)) {
    estimated$params$beta = NULL
    MSE_beta = 0
  }
  
  if (ncol(estimated$params$Lambda) >= ncol(true$Lambda)) {
    combos= combn(ncol(estimated$params$Lambda), ncol(true$Lambda))
    num_combos = ncol(combos)
    tmp = numeric(num_combos)
    for (i in 1:num_combos) {
      current_idx = combos[,i]
      A_sub = estimated$params$Lambda[, current_idx]
      tmp[i] = procrustes(A_sub, true$Lambda, symmetric = TRUE)$ss
    }
    best_idx_in_combos = which.min(tmp)
    best_cols = combos[, best_idx_in_combos]
    errors2 = evaluate_Ex4(true$Lambda, estimated$params$Lambda[, best_cols])
  } else {
    errors2 = evaluate_Ex4(true$Lambda, estimated$params$Lambda)
  }
  
  
  if (ncol(estimated$params$beta) == ncol(true$beta)) {
    MSE_beta = MSE(estimated$params$beta, true$beta)
  } else {
    MSE_beta = 0
  }
  
  mat1 = cbind(estimated$params$B, estimated$params$beta, estimated$params$Lambda)
  mat2 = cbind(true$B, true$beta, true$Lambda)
  
  
  
  if (ncol(mat1) == ncol(mat2)) {
    errors1 = evaluate_Ex4(mat2, mat1)
  } else if (ncol(mat1) > ncol(mat2) ) {
    if (ncol(estimated$params$beta) > ncol(true$beta) && (ncol(estimated$params$Lambda) == ncol(true$Lambda))) {
      errors1 = evaluate_Ex4(cbind(estimated$params$B, estimated$params$beta[,1:ncol(true$beta)], estimated$params$Lambda), mat2)
    } else if (ncol(estimated$params$beta) == ncol(true$beta) && (ncol(estimated$params$Lambda) > ncol(true$Lambda))) {
      errors1 = evaluate_Ex4(cbind(estimated$params$B, estimated$params$beta, estimated$params$Lambda[,best_cols]), mat2)
    } else if (ncol(estimated$params$beta) > ncol(true$beta) && (ncol(estimated$params$Lambda) > ncol(true$Lambda))) {
      errors1 = evaluate_Ex4(cbind(estimated$params$B, estimated$params$beta[,1:ncol(true$beta)], estimated$params$Lambda[,best_cols]), mat2)
    }
  } else {
    errors1 = evaluate_Ex4(mat2, mat1)
  }
  
  c(errors1, errors2, MSE_pai, MSE_B, MSE_beta, MSE_emat, MSE_ll)
}



evaluate_FACUB3 = function(true, estimated) {
  MSE_pai = MSE(estimated$params$pai, true$pai)
  MSE_B = MSE(estimated$params$B, true$B)
  MSE_emat = MSE(estimated$e.mat, true$csi)
  MSE_ll = MSE(estimated$ll, true$csi)
  if (is.null(true$beta)) {
    estimated$params$beta = NULL
    MSE_beta = 0
  }
  
  if (ncol(estimated$params$Lambda) >= ncol(true$Lambda)) {
    combos= combn(ncol(estimated$params$Lambda), ncol(true$Lambda))
    num_combos = ncol(combos)
    tmp = numeric(num_combos)
    for (i in 1:num_combos) {
      current_idx = combos[,i]
      A_sub = estimated$params$Lambda[, current_idx]
      tmp[i] = procrustes(A_sub, true$Lambda, symmetric = TRUE)$ss
    }
    best_idx_in_combos = which.min(tmp)
    best_cols = combos[, best_idx_in_combos]
    errors2 = evaluate_Ex4(true$Lambda, estimated$params$Lambda[, best_cols])
  } else {
    errors2 = rep(0,4)
  }
  
  
  if (ncol(estimated$params$beta) == ncol(true$beta)) {
    MSE_beta = MSE(estimated$params$beta, true$beta)
  } else if (ncol(estimated$params$beta) > ncol(true$beta)) {
    MSE_beta = MSE(estimated$params$beta[,1:ncol(true$beta)], true$beta)
  } else {
    MSE_beta = MSE(estimated$params$beta, true$beta[,1:ncol(estimated$params$beta)])
  }
  
  mat1 = cbind(estimated$params$B, estimated$params$beta, estimated$params$Lambda)
  mat2 = cbind(true$B, true$beta, true$Lambda)
  
  
  
  if (ncol(mat1) == ncol(mat2)) {
    errors1 = evaluate_Ex4(mat2, mat1)
  } else if (ncol(mat1) > ncol(mat2) ) {
    if (ncol(estimated$params$beta) > ncol(true$beta) && (ncol(estimated$params$Lambda) == ncol(true$Lambda))) {
      errors1 = evaluate_Ex4(cbind(estimated$params$B, estimated$params$beta[,1:ncol(true$beta)], estimated$params$Lambda), mat2)
    } else if (ncol(estimated$params$beta) == ncol(true$beta) && (ncol(estimated$params$Lambda) > ncol(true$Lambda))) {
      errors1 = evaluate_Ex4(cbind(estimated$params$B, estimated$params$beta, estimated$params$Lambda[,best_cols]), mat2)
    } else if (ncol(estimated$params$beta) > ncol(true$beta) && (ncol(estimated$params$Lambda) > ncol(true$Lambda))) {
      errors1 = evaluate_Ex4(cbind(estimated$params$B, estimated$params$beta[,1:ncol(true$beta)], estimated$params$Lambda[,best_cols]), mat2)
    }
  } else {
    errors1 = c(procrustes(mat1, mat2, symmetric=TRUE)$ss, rep(0,3))
  }
  
  c(errors1, errors2, MSE_pai, MSE_B, MSE_beta, MSE_emat, MSE_ll)
}


evaluate_MCUB = function(true, estimated) {
  MSE_pai = MSE(estimated$pai, true$pai)
  MSE_B = MSE(estimated$B, true$B)
  if (is.null(true$beta)) {
    estimated$params$beta = NULL
    MSE_beta = 0
  } else {
    MSE_beta = MSE(estimated$beta, true$beta)
  }
  MSE_emat = MSE_ll = MSE(estimated$csi, true$csi)
  mat1 = cbind(true$B, true$beta, true$Lambda)
  mat2 = cbind(estimated$B, estimated$beta)
  c(evaluate_Ex4(mat1, mat2), MSE_pai, MSE_B, MSE_beta, MSE_emat, MSE_ll)
}

evaluate_MCUBo = function(true, estimated) {
  p = nrow(true$Lambda)
  d = ncol(true$Lambda)
  MSE_pai = MSE(estimated$pai, true$pai)
  MSE_B = MSE(estimated$B, true$B)
  if (is.null(true$beta)) {
    estimated$params$beta = NULL
    MSE_beta = 0
  } else {
    MSE_beta = MSE(estimated$beta[, 1:ncol(true$beta)], true$beta)
  }
  MSE_emat = MSE_ll = MSE(estimated$csi, true$csi)
  mat1 = cbind(true$B, true$beta, true$Lambda)
  mat2 = cbind(estimated$B, estimated$beta)
  c(evaluate_Ex4(mat1, mat2), MSE_pai, MSE_B, MSE_beta, MSE_emat, MSE_ll)
}


generate_data_lavaan = function(n, p, d) {
  Phi = matrix(c(1, .2, .5, .2, 1, .8, .5, .8, 1), 3, 3)
  choice = seq(0.3, 0.9, 0.1)
  Lambda = matrix(c(sample(choice, 17, replace = TRUE), rep(0, p-1), sample(choice, 17, replace = TRUE), rep(0, p-1), sample(choice, 18, replace = TRUE)), p, d)
  f = rmvn(n, rep(0, d), Phi)
  #Lambda = matrix(runif(p*d, 0, 1), p, d)
  Theta = diag(p) - diag(diag(Lambda %*% Phi %*% t(Lambda)))
  delta = rmvn(n, rep(0, p), Theta)
  #f = matrix(NA, n, d)
  #for (i in 1:n) {
  #  f[i, ] = rnorm(d, 0, 1)
  #}
  xstar = f %*% t(Lambda) + delta
  ind1 = which(xstar < -1.2)
  ind2 = which(xstar >= -1.2 & xstar < 0)
  ind3 = which(xstar >= 0 & xstar < 1.2)
  ind4 = which(xstar >= 1.2)
  x = matrix(0, n, p)
  x[ind1] = 1
  x[ind2] = 2
  x[ind3] = 3
  x[ind4] = 4
  return(list(R = x, f = f, Lambda = Lambda))
}

generate_data_lavaan2 = function(n, p, d) {
  Phi = matrix(c(1, .2, .5, .2, 1, .8, .5, .8, 1), 3, 3)
  choice = seq(0.3, 0.9, 0.1)
  Lambda = matrix(c(sample(choice, 17, replace = TRUE), rep(0, p-1), sample(choice, 17, replace = TRUE), rep(0, p-1), sample(choice, 18, replace = TRUE)), p, d)
  f = rmvn(n, rep(0, d), Phi)
  #Lambda = matrix(runif(p*d, 0, 1), p, d)
  Theta = diag(p) - diag(diag(Lambda %*% Phi %*% t(Lambda)))
  delta = rmvn(n, rep(0, p), Theta)
  #f = matrix(NA, n, d)
  #for (i in 1:n) {
  #  f[i, ] = rnorm(d, 0, 1)
  #}
  xstar = f %*% t(Lambda) + delta
  ind1 = which(xstar < -3)
  ind2 = which(xstar >= -3 & xstar < -1)
  ind3 = which(xstar >= -1 & xstar < 1)
  ind4 = which(xstar >= 1 & xstar < 3)
  ind5 = which(xstar >= 3)
  x = matrix(0, n, p)
  x[ind1] = 1
  x[ind2] = 2
  x[ind3] = 3
  x[ind4] = 4
  x[ind5] = 5
  return(list(R = x, f = f, Lambda = Lambda))
}


FAVA_old = function (R, V = NULL, m, n.factors, init = NULL, maxit = 200, trace = FALSE) 
{
  n = nrow(R)
  p = ncol(R)
  if (is.null(V)) {
    q = 1
    X = matrix(0, n, q)
  } else {
    X = V
    q = ncol(X)
  }
  out.list = list()
  R_sc = scale(log2(1 + R), center = TRUE, scale = TRUE)
  re = svd(R_sc, n.factors, n.factors)
  if (n.factors == 1) {
    mu = new.mu = re$u * (re$d[1])
  } else {
    mu = new.mu = re$u %*% diag(re$d[1:n.factors])
  }
  Lambda = new.Lambda = re$v
  #mu = new.mu = res3$lvs$mu
  #Lambda = new.Lambda = matrix(c(.9, .8 , .7, .5, rep(0, 5), .6, .7, .8), p, d)
  B = new.B = rep(1, p)
  beta = new.beta = matrix(1, p, q)
  factor_coefs = cbind(B, beta, Lambda)
  sigma = new.sigma = matrix(0.01, n, n.factors)
  #sigma = new.sigma = res3$lvs$sigma
  pai = new.pai = rep(0.5, p)
  alpha = new.alpha = matrix(pai, n, p, byrow = TRUE)
  if (!is.null(init)) {
    mu=new.mu=init$mu
    sigma=new.sigma=init$sigma
    B=new.B=init$B
    beta=new.beta=init$beta
    Lambda=new.Lambda=init$Lambda
    pai=new.pai=init$pai
    alpha=new.alpha=init$alpha
  }
  cur.VLB = -1e+06
  iter = 1
  ratio = 10
  diff = 1e+05
  eps = 1e-04
  max.iter = 100
  m.cur.logfunc = -1e+06
  b.cur.logfunc = -1e+06
  d.cur.logfunc = -1e+06
  s.cur.logfunc = -1e+06
  while ((diff > eps * (abs(cur.VLB) + eps)) && iter <= max.iter) {
    if (trace) 
      cat("Iteration: ", iter, "\n")
    VLB = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      fun1 = function(i) {
        -0.5 * (sum(new.sigma[i, ]) + sum(new.mu[i, ]^2) - 
                  sum(log(new.sigma[i, ])))
      }
      fun2 = function(i) {
        new.pai = new.pai + 1e-08
        new.alpha = new.alpha + 1e-08
        alpha.mat = (1 - new.alpha[i, ]) * log((1 - new.pai)/(1 - 
                                                                new.alpha[i, ])) + new.alpha[i, ] * log(new.pai/new.alpha[i, 
                                                                ])
        sum(alpha.mat)
      }
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      fun3 = new.alpha * (log(choose(m - 1, R - 1)) + (R - 
                                                         1) * ll - (m - 1) * log(1 + exp(e.mat))) - (1 - 
                                                                                                       new.alpha) * log(m)
      y = sum(sapply(1:n, fun1)) + sum(sapply(1:n, fun2)) + 
        sum(fun3)
      return(y)
    }
    log_base = function(model.coefs, va.mu, va.sigma, pai, 
                        alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      y1 = new.alpha * (log(choose(m - 1, R - 1)) + (R - 
                                                       1) * ll - (m - 1) * log(1 + exp(ll))) - (1 - 
                                                                                                  new.alpha) * log(m)
      return(sum(y1))
    }
    log_base2 = function(model.coefs, va.mu, va.sigma, pai, 
                         alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      y1 = new.alpha * (log(choose(m - 1, R - 1)) + (R - 
                                                       1) * ll - (m - 1) * log(1 + exp(e.mat))) - (1 - 
                                                                                                     new.alpha) * log(m)
      return(sum(y1))
    }
    new.pai = round(apply(new.alpha, 2, mean), 6)
    model.coefs_f = function(model.coefs, va.mu, va.sigma, 
                             pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      y1 = new.alpha * (R - 1) * ll
      y2 = -new.alpha * (m - 1) * log(1 + exp(e.mat))
      y = sum(y1) + sum(y2)
      return(y)
    }
    model.coefs_gr_f = function(model.coefs, va.mu, va.sigma, 
                                pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      b1 = new.alpha * (R - 1) - (m - 1) * new.alpha * 
        exp(e.mat)/(1 + exp(e.mat))
      b2 = NULL
      for (w in 1:n.factors) {
        b2 = c(b2, new.alpha * (sweep(R - 1, 1, new.mu[, 
                                                       w], "*") - (m - 1) * sweep((new.sigma[, w] %*% 
                                                                                     t(new.Lambda[, w])), 1, new.mu[, w], "+") * 
                                  exp(e.mat)/(1 + exp(e.mat))))
      }
      b2 = matrix(b2, n, p * n.factors)
      b3 = NULL
      for (w in 1:q) {
        b3 = c(b3, new.alpha * (sweep(R - 1 - (m - 1) * 
                                        exp(e.mat)/(1 + exp(e.mat)), 1, X[, w], "*")))
      }
      b3 = matrix(b3, n, p * q)
      return(c(colSums(b2), colSums(b1), colSums(b3)))
    }
    qq = try(optim(c(Lambda, B, beta), va.mu = c(new.mu), 
                   va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha, 
                   method = "BFGS", fn = model.coefs_f, gr = model.coefs_gr_f, 
                   control = list(trace = 0, fnscale = -1, maxit = maxit)), 
             silent = TRUE)
    if ("try-error" %in% class(qq)) {
      new.Lambda = Lambda
      new.B = B
      new.beta = beta
    } else {
      if (iter > 1 && b.cur.logfunc > qq$value) {
        if (trace) {
          cat("Optimization of model coefs dit not improve on iteration step ", 
              iter, "\n")
        }
        new.Lambda = Lambda
        new.B = B
        new.beta = beta
      } else {
        if (trace) {
          cat("Model parameters updated", "\n")
        }
        new.Lambda = matrix(qq$par[1:(p * n.factors)], 
                            p, n.factors)
        qq$par = qq$par[-(1:(p * n.factors))]
        new.B = qq$par[1:p]
        qq$par = qq$par[-(1:p)]
        new.beta = matrix(qq$par[1:(p * q)], p, q)
        if (qq$convergence != 0) {
          if (trace) {
            cat("Optimization of model coefs did not converge on iteration step ", 
                iter, "\n")
          }
        }
      }
    }
    delta.alpha.required = 0.001
    p.iter = 1
    p.max.iter = 10
    delta.alpha = abs(alpha)
    while (!all(delta.alpha < delta.alpha.required) & (p.iter < 
                                                       p.max.iter)) {
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      new.alpha = matrix(0, n, p)
      for (i in 1:n) {
        new.alpha[i, ] = plogis(log(new.pai/(1 - new.pai)) + 
                                  log(m) + log(choose(m - 1, R[i, ] - 1)) + (R[i, 
                                  ] - 1) * ll[i, ] - (m - 1) * log(1 + exp(e.mat[i, 
                                  ])))
      }
      va.sigma_f = function(model.coefs, va.mu, va.sigma, 
                            pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% 
          t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        fun1 = function(i) {
          -0.5 * (sum(new.sigma[i, ]) - sum(log(new.sigma[i, 
          ])))
        }
        y2 = -(m - 1) * new.alpha * log(1 + exp(e.mat))
        y = sum(sapply(1:n, fun1)) + sum(y2)
        return(y)
      }
      va.sigma_gr_f = function(model.coefs, va.mu, va.sigma, 
                               pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% 
          t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        beta2 = new.Lambda^2
        mu.mat = -(m - 1) * new.alpha * exp(e.mat)/(1 + 
                                                      exp(e.mat))
        grad.sigma = matrix(NA, n, n.factors)
        for (i in 1:n) {
          grad.sigma[i, ] = diag(-0.5 * (diag(rep(1, 
                                                  n.factors)) - (diag(new.sigma[i, ]))^(-1)) + 
                                   diag(apply(mu.mat[i, ] * beta2, 2, sum)))
        }
        return(c(grad.sigma))
      }
      qq = try(constrOptim(c(sigma), method = "BFGS", f = va.sigma_f, 
                           gr = va.sigma_gr_f, model.coefs = c(new.Lambda, 
                                                               new.B, new.beta), va.mu = c(new.mu), pai = new.pai, 
                           alpha = new.alpha, ui = diag(1, n * n.factors, 
                                                        n * n.factors), ci = rep(1e-08, n * n.factors), 
                           control = list(trace = 0, fnscale = -1, maxit = maxit, 
                                          reltol = 1e-06)), silent = TRUE)
      if ("try-error" %in% class(qq)) {
        new.sigma = sigma
      } else {
        if (iter > 1 && s.cur.logfunc > qq$value) {
          if (trace) {
            cat("Optimization of sigma did not improve on iteration step ", 
                iter, "\n")
          }
          new.sigma = sigma
        } else {
          if (trace) {
            cat("Variational parameters sigma updated", 
                "\n")
          }
          new.sigma = matrix(qq$par, n, n.factors)
          if (qq$convergence != 0) {
            if (trace) {
              cat("Optimization of sigma did not converge on iteration step ", 
                  iter, "\n")
            }
          }
        }
      }
      #if (!is.null(init)) new.sigma=init$sigma
      mu_f = function(model.coefs, va.mu, va.sigma, pai, 
                      alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        fun1 = function(i) {
          sum(new.mu[i, ]^2)
        }
        y2 = new.alpha * (R - 1) * new.mu %*% t(new.Lambda)
        y3 = -(m - 1) * new.alpha * log(1 + exp(e.mat))
        y = sum(sapply(1:n, fun1)) + sum(y2) + sum(y3)
        return(y)
      }
      mu_grad_f = function(model.coefs, va.mu, va.sigma, 
                           pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% 
          t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        grad.m = NULL
        sum1 = new.alpha * ((R - 1) - (m - 1) * exp(e.mat)/(1 + 
                                                              exp(e.mat)))
        for (w in 1:n.factors) {
          grad.m = c(grad.m, rowSums(sweep(sum1, 2, new.Lambda[, 
                                                               w], "*")) - new.mu[, w])
        }
        return(c(grad.m))
      }
      qq = try(optim(c(mu), method = "BFGS", fn = mu_f, 
                     gr = mu_grad_f, model.coefs = c(new.Lambda, new.B, 
                                                     new.beta), va.sigma = c(new.sigma), pai = new.pai, 
                     alpha = new.alpha, control = list(trace = 0, 
                                                       fnscale = -1, maxit = maxit, reltol = 1e-06)), 
               silent = TRUE)
      if ("try-error" %in% class(qq)) {
        new.mu = mu
      } else {
        if (iter > 1 && m.cur.logfunc > qq$value) {
          if (trace) {
            cat("Optimization of mu did not improve on iteration step ", 
                iter, "\n")
          }
          new.mu = mu
        } else {
          if (trace) {
            cat("Variational parameters mu updated", 
                "\n")
          }
          new.mu = matrix(qq$par, n, n.factors)
          if (qq$convergence != 0) {
            if (trace) {
              cat("Optimization of mu did not converge on iteration step ", 
                  iter, "\n")
            }
          }
        }
      }
      #if (!is.null(init)) new.mu=init$mu
      q_m = list(value = mu_f(c(new.Lambda, new.B, new.beta), 
                              va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, 
                              alpha = new.alpha))
      m.new.logfunc = q_m$value
      m.cur.logfunc = m.new.logfunc
      q_s = list(value = va.sigma_f(c(new.Lambda, new.B, 
                                      new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), 
                                    pai = new.pai, alpha = new.alpha))
      s.new.logfunc = q_s$value
      s.cur.logfunc = s.new.logfunc
      delta.alpha = abs(new.alpha - alpha)
      sigma = new.sigma
      mu = new.mu
      alpha = new.alpha
      pai = new.pai
      p.iter = p.iter + 1
    }
    q_b = list(value = model.coefs_f(c(new.Lambda, new.B, 
                                       new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), 
                                     pai = new.pai, alpha = new.alpha))
    b.new.logfunc = q_b$value
    b.cur.logfunc = b.new.logfunc
    qq = list(value = VLB(c(new.Lambda, new.B, new.beta), 
                          va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, 
                          alpha = new.alpha))
    new.VLB = qq$value
    diff = abs(new.VLB - cur.VLB)
    ratio = abs(new.VLB/cur.VLB)
    if (trace) 
      cat("New VLB: ", new.VLB, "cur VLB: ", cur.VLB, "Ratio of VLB: ", 
          ratio, ". Difference in VLB: ", diff, "\n")
    #if (trace) cat("\nThe current evaluate results: ", evaluate_Ex4(promax(new.Lambda)$loadings, promax(data$Lambda)$loadings))
    cur.VLB = new.VLB
    qq = list(value = log_base(c(new.Lambda, new.B, new.beta), 
                               va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, 
                               alpha = new.alpha))
    cur.log = qq$value
    EE = log_base2(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), 
                   va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha)
    EN2 = 2 * cur.log - 2 * EE
    pai = new.pai
    alpha = new.alpha
    B = new.B
    Lambda = new.Lambda
    beta = new.beta
    mu = new.mu
    sigma = new.sigma
    iter = iter + 1
  }
  if (iter > 99) {
    print("FACUB not converging!")
  }
  ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
    new.mu %*% t(new.Lambda)
  e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
  if (is.null(V)) beta = NULL
  out.list$VLB = cur.VLB
  out.list$EE = EE
  out.list$lob = cur.log
  out.list$EN2 = EN2
  out.list$iter = iter - 1
  out.list$lvs$alpha = alpha
  out.list$lvs$mu = mu
  out.list$lvs$sigma = sigma
  out.list$params$pai = pai
  out.list$params$B = B
  out.list$params$beta = beta
  out.list$params$Lambda = Lambda
  out.list$ll = plogis(ll)
  out.list$e.mat = plogis(e.mat)
  return(out.list)
}



FAVA_old2 = function (R, V = NULL, m, n.factors, init = NULL, maxit = 200, trace = FALSE) 
{
  n = nrow(R)
  p = ncol(R)
  if (is.null(V)) {
    q = 1
    X = matrix(0, n, q)
  } else {
    X = V
    q = ncol(X)
  }
  out.list = list()
  R_sc = scale(log2(1 + R), center = TRUE, scale = TRUE)
  re = svd(R_sc, n.factors, n.factors)
  if (n.factors == 1) {
    mu = new.mu = re$u * (re$d[1])
  } else {
    mu = new.mu = re$u %*% diag(re$d[1:n.factors])
  }
  Lambda = new.Lambda = re$v
  #mu = new.mu = res3$lvs$mu
  #Lambda = new.Lambda = matrix(c(.9, .8 , .7, .5, rep(0, 5), .6, .7, .8), p, d)
  B = new.B = rep(1, p)
  beta = new.beta = matrix(1, p, q)
  factor_coefs = cbind(B, beta, Lambda)
  sigma = new.sigma = matrix(0.01, n, n.factors)
  #sigma = new.sigma = res3$lvs$sigma
  pai = new.pai = rep(0.5, p)
  alpha = new.alpha = matrix(pai, n, p, byrow = TRUE)
  if (!is.null(init)) {
    mu=new.mu=init$mu
    sigma=new.sigma=init$sigma
    B=new.B=init$B
    beta=new.beta=init$beta
    Lambda=new.Lambda=init$Lambda
    pai=new.pai=init$pai
    alpha=new.alpha=init$alpha
  }
  cur.VLB = -1e+06
  iter = 1
  ratio = 10
  diff = 1e+05
  eps = 1e-04
  max.iter = 100
  m.cur.logfunc = -1e+06
  b.cur.logfunc = -1e+06
  d.cur.logfunc = -1e+06
  s.cur.logfunc = -1e+06
  while ((diff > eps * (abs(cur.VLB) + eps)) && iter <= max.iter) {
    if (trace) 
      cat("Iteration: ", iter, "\n")
    VLB = function(model.coefs, va.mu, va.sigma, pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      fun1 = function(i) {
        -0.5 * (sum(new.sigma[i, ]) + sum(new.mu[i, ]^2) - 
                  sum(log(new.sigma[i, ])))
      }
      fun2 = function(i) {
        new.pai = new.pai + 1e-08
        new.alpha = new.alpha + 1e-08
        alpha.mat = (1 - new.alpha[i, ]) * log((1 - new.pai)/(1 - 
                                                                new.alpha[i, ])) + new.alpha[i, ] * log(new.pai/new.alpha[i, 
                                                                ])
        sum(alpha.mat)
      }
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      fun3 = new.alpha * (log(choose(m - 1, R - 1)) + (R - 
                                                         1) * ll - (m - 1) * log(1 + exp(e.mat))) - (1 - 
                                                                                                       new.alpha) * log(m)
      y = sum(sapply(1:n, fun1)) + sum(sapply(1:n, fun2)) + 
        sum(fun3)
      return(y)
    }
    log_base = function(model.coefs, va.mu, va.sigma, pai, 
                        alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      y1 = new.alpha * (log(choose(m - 1, R - 1)) + (R - 
                                                       1) * ll - (m - 1) * log(1 + exp(ll))) - (1 - 
                                                                                                  new.alpha) * log(m)
      return(sum(y1))
    }
    log_base2 = function(model.coefs, va.mu, va.sigma, pai, 
                         alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      y1 = new.alpha * (log(choose(m - 1, R - 1)) + (R - 
                                                       1) * ll - (m - 1) * log(1 + exp(e.mat))) - (1 - 
                                                                                                     new.alpha) * log(m)
      return(sum(y1))
    }
    new.pai = round(apply(new.alpha, 2, mean), 6)
    model.coefs_f = function(model.coefs, va.mu, va.sigma, 
                             pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      y1 = new.alpha * (R - 1) * ll
      y2 = -new.alpha * (m - 1) * log(1 + exp(e.mat))
      y = sum(y1) + sum(y2)
      return(y)
    }
    model.coefs_gr_f = function(model.coefs, va.mu, va.sigma, 
                                pai, alpha) {
      new.mu = matrix(va.mu, n, n.factors)
      new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                          p, n.factors)
      model.coefs = model.coefs[-(1:(p * n.factors))]
      new.B = model.coefs[1:p]
      model.coefs = model.coefs[-(1:p)]
      new.beta = matrix(model.coefs[1:(p * q)], p, q)
      new.sigma = matrix(va.sigma, n, n.factors)
      new.pai = pai
      new.alpha = alpha
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      b1 = new.alpha * (R - 1) - (m - 1) * new.alpha * 
        exp(e.mat)/(1 + exp(e.mat))
      b2 = NULL
      for (w in 1:n.factors) {
        b2 = c(b2, new.alpha * (sweep(R - 1, 1, new.mu[, 
                                                       w], "*") - (m - 1) * sweep((new.sigma[, w] %*% 
                                                                                     t(new.Lambda[, w])), 1, new.mu[, w], "+") * 
                                  plogis(e.mat)))
      }
      b2 = matrix(b2, n, p * n.factors)
      b3 = NULL
      for (w in 1:q) {
        b3 = c(b3, new.alpha * (sweep(R - 1 - (m - 1) * plogis(e.mat), 1, X[, w], "*")))
      }
      b3 = matrix(b3, n, p * q)
      return(c(colSums(b2), colSums(b1), colSums(b3)))
    }
    qq = try(optim(c(Lambda, B, beta), va.mu = c(new.mu), 
                   va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha, 
                   method = "BFGS", fn = model.coefs_f, gr = model.coefs_gr_f, 
                   control = list(trace = 0, fnscale = -1, maxit = maxit)), 
             silent = TRUE)
    if ("try-error" %in% class(qq)) {
      new.Lambda = Lambda
      new.B = B
      new.beta = beta
    } else {
      if (iter > 1 && b.cur.logfunc > qq$value) {
        if (trace) {
          cat("Optimization of model coefs dit not improve on iteration step ", 
              iter, "\n")
        }
        new.Lambda = Lambda
        new.B = B
        new.beta = beta
      } else {
        if (trace) {
          cat("Model parameters updated", "\n")
        }
        new.Lambda = matrix(qq$par[1:(p * n.factors)], 
                            p, n.factors)
        qq$par = qq$par[-(1:(p * n.factors))]
        new.B = qq$par[1:p]
        qq$par = qq$par[-(1:p)]
        new.beta = matrix(qq$par[1:(p * q)], p, q)
        if (qq$convergence != 0) {
          if (trace) {
            cat("Optimization of model coefs did not converge on iteration step ", 
                iter, "\n")
          }
        }
      }
    }
    delta.alpha.required = 0.001
    p.iter = 1
    p.max.iter = 10
    delta.alpha = abs(alpha)
    while (!all(delta.alpha < delta.alpha.required) & (p.iter < 
                                                       p.max.iter)) {
      ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
        new.mu %*% t(new.Lambda)
      e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
      new.alpha = matrix(0, n, p)
      # update variational parameters: alpha
      for (i in 1:n) {
        new.alpha[i, ] = plogis(log(new.pai/(1 - new.pai)) + 
                                  log(m) + log(choose(m - 1, R[i, ] - 1)) + (R[i, 
                                  ] - 1) * ll[i, ] - (m - 1) * log(1 + exp(e.mat[i, 
                                  ])))
      }
      # update variational parameters: sigma
      va.sigma_f = function(model.coefs, va.mu, va.sigma, 
                            pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% 
          t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        fun1 = function(i) {
          -0.5 * (sum(new.sigma[i, ]) - sum(log(new.sigma[i, 
          ])))
        }
        y2 = -(m - 1) * new.alpha * log(1 + exp(e.mat))
        y = sum(sapply(1:n, fun1)) + sum(y2)
        return(y)
      }
      
      va.sigma_gr_f = function(model.coefs, va.mu, va.sigma, 
                               pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% 
          t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        beta2 = new.Lambda^2
        mu.mat = -(m - 1) * new.alpha * exp(e.mat)/(1 + 
                                                      exp(e.mat))
        grad.sigma = matrix(NA, n, n.factors)
        for (i in 1:n) {
          grad.sigma[i, ] = diag(-0.5 * (diag(rep(1, 
                                                  n.factors)) - diag(1/new.sigma[i, ], n.factors)) + 
                                   diag(apply(mu.mat[i, ] * beta2, 2, sum)))
        }
        return(c(grad.sigma))
      }
      
      va.sigma_fi = function(model.coefs, va.mu, va.sigma, 
                             pai, alpha, ii) {
        new.mu = va.mu # n.factors dimension
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = va.sigma # n.factors dimension
        new.pai = pai
        new.alpha = alpha # p dimension
        ll = new.B + X[ii, ] %*% t(new.beta) + new.mu %*% t(new.Lambda) # p dimension
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        y1 = -0.5 * (sum(new.sigma) - sum(log(new.sigma)))
        y2 = -(m - 1) * new.alpha * log(1 + exp(e.mat))
        y = y1 + sum(y2)
        return(y)
      }
      
      va.sigma_gr_fi = function(model.coefs, va.mu, va.sigma, 
                                pai, alpha, ii) {
        new.mu = va.mu # n.factors dimension
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = va.sigma # n.factors dimension
        new.pai = pai
        new.alpha = alpha # p dimension
        ll = c(new.B + X[ii, ] %*% t(new.beta) + new.mu %*% t(new.Lambda)) # p dimension
        e.mat = c(ll + 0.5 * (new.sigma) %*% t(new.Lambda^2))
        beta2 = new.Lambda^2
        mu.vec = -(m - 1) * new.alpha * plogis(e.mat)
        grad.sigma = diag(-0.5 * (diag(rep(1, n.factors), n.factors) - diag(new.sigma^(-1), n.factors)) + 
                            diag(colSums(mu.vec * beta2), n.factors))
        
        return(c(grad.sigma))
      }
      for (i in 1:n) {
        qq = try(constrOptim(sigma[i, ], method = "BFGS", f = va.sigma_fi,
                             gr = va.sigma_gr_fi,
                             model.coefs = c(new.Lambda,new.B, new.beta),
                             va.mu = new.mu[i, ],
                             pai = new.pai,
                             alpha = new.alpha[i, ],
                             ii = i,
                             ui = diag(1, n.factors, n.factors),
                             ci = rep(1e-08, n.factors),
                             control = list(trace = 0,
                                            fnscale = -1,
                                            maxit = maxit,
                                            reltol = 1e-06)), silent = TRUE)
        if ("try-error" %in% class(qq)) {
          new.sigma[i, ] = sigma[i, ]
        } else {
          if (iter > 1 && s.cur.logfunc > qq$value) {
            if (trace) {
              cat("Optimization of sigma did not improve on iteration step ",
                  iter, "\n")
            }
            new.sigma[i, ] = sigma[i, ]
          } else {
            if (trace) {
              cat("Variational parameters sigma updated",
                  "\n")
            }
            new.sigma[i, ] = qq$par
            if (qq$convergence != 0) {
              if (trace) {
                cat("Optimization of sigma did not converge on iteration step ",
                    iter, "\n")
              }
            }
          }
        }
      }
      # qq = try(constrOptim(c(sigma), method = "BFGS", 
      #                      f = va.sigma_f, 
      #                      gr = va.sigma_gr_f, 
      #                      model.coefs = c(new.Lambda, new.B, new.beta), 
      #                      va.mu = c(new.mu), 
      #                      pai = new.pai, 
      #                      alpha = new.alpha, 
      #                      ui = diag(1, n * n.factors, n * n.factors), 
      #                      ci = rep(1e-08, n * n.factors), 
      #                      control = list(trace = 0, fnscale = -1, maxit = maxit, 
      #                                     reltol = 1e-06)), silent = TRUE)
      # if ("try-error" %in% class(qq)) {
      #   new.sigma = sigma
      # } else {
      #   if (iter > 1 && s.cur.logfunc > qq$value) {
      #     if (trace) {
      #       cat("Optimization of sigma did not improve on iteration step ", 
      #           iter, "\n")
      #     }
      #     new.sigma = sigma
      #   } else {
      #     if (trace) {
      #       cat("Variational parameters sigma updated", 
      #           "\n")
      #     }
      #     new.sigma = matrix(qq$par, n, n.factors)
      #     if (qq$convergence != 0) {
      #       if (trace) {
      #         cat("Optimization of sigma did not converge on iteration step ", 
      #             iter, "\n")
      #       }
      #     }
      #   }
      # }
      #if (!is.null(init)) new.sigma=init$sigma
      mu_f = function(model.coefs, va.mu, va.sigma, pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        fun1 = function(i) {
          sum(new.mu[i, ]^2)
        }
        y2 = new.alpha * (R - 1) * new.mu %*% t(new.Lambda)
        y3 = -(m - 1) * new.alpha * log(1 + exp(e.mat))
        y = -0.5*sum(sapply(1:n, fun1)) + sum(y2) + sum(y3)
        return(y)
      }
      mu_grad_f = function(model.coefs, va.mu, va.sigma, 
                           pai, alpha) {
        new.mu = matrix(va.mu, n, n.factors)
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = matrix(va.sigma, n, n.factors)
        new.pai = pai
        new.alpha = alpha
        ll = matrix(new.B, n, p, byrow = TRUE) + X %*% 
          t(new.beta) + new.mu %*% t(new.Lambda)
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        grad.m = NULL
        sum1 = new.alpha * ((R - 1) - (m - 1) * plogis(e.mat))
        for (w in 1:n.factors) {
          grad.m = c(grad.m, rowSums(sweep(sum1, 2, new.Lambda[, w], "*")) - new.mu[, w])
        }
        return(c(grad.m))
      }
      mu_fi = function(model.coefs, va.mu, va.sigma, pai, alpha, ii) {
        new.mu = va.mu # n.factors dimension
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p, 
                          q)
        new.sigma = va.sigma # n.factors dimension
        new.pai = pai
        new.alpha = alpha # p dimension
        ll = new.B + X[ii, ] %*% t(new.beta) + new.mu %*% t(new.Lambda) # p dimension
        e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
        y1 = -0.5*sum(new.mu^2)
        y2 = new.alpha * (R[ii,] - 1) * new.mu %*% t(new.Lambda)
        y3 = -(m - 1) * new.alpha * log(1 + exp(e.mat))
        y = y1 + sum(y2) + sum(y3)
        return(y)
      }
      mu_grad_fi = function(model.coefs, va.mu, va.sigma, pai, alpha, ii) {
        new.mu = va.mu # n.factors dimension
        new.Lambda = matrix(model.coefs[1:(p * n.factors)], 
                            p, n.factors)
        model.coefs = model.coefs[-(1:(p * n.factors))]
        new.B = model.coefs[1:p]
        model.coefs = model.coefs[-(1:p)]
        new.beta = matrix(model.coefs[1:(p * q)], p,q)
        new.sigma = va.sigma # n.factors dimension
        new.pai = pai
        new.alpha = alpha # p dimension
        ll = c(new.B + X[ii, ] %*% t(new.beta) + new.mu %*% t(new.Lambda)) # p dimension
        e.mat = c(ll + 0.5 * (new.sigma) %*% t(new.Lambda^2))
        sum1 = new.alpha * ((R[ii, ] - 1) - (m - 1) * plogis(e.mat))
        return(-new.mu + as.vector(crossprod(sum1, new.Lambda)))
      }
      for (i in 1:n) {
        qq = try(optim(mu[i, ], method = "BFGS", fn = mu_fi,
                       gr = mu_grad_fi,
                       model.coefs = c(new.Lambda, new.B, new.beta),
                       va.sigma = new.sigma[i, ],
                       pai = new.pai,
                       alpha = new.alpha[i, ],
                       ii = i,
                       control = list(trace = 0,
                                      fnscale = -1,
                                      maxit = maxit,
                                      reltol = 1e-06)), silent = TRUE)
        if ("try-error" %in% class(qq)) {
          new.mu[i, ] = mu[i, ]
        } else {
          if (iter > 1 && m.cur.logfunc > qq$value) {
            if (trace) {
              cat("Optimization of mu did not improve on iteration step ",
                  iter, "\n")
            }
            new.mu[i, ] = mu[i, ]
          } else {
            if (trace) {
              cat("Variational parameters mu updated",
                  "\n")
            }
            new.mu[i, ] = qq$par
            if (qq$convergence != 0) {
              if (trace) {
                cat("Optimization of mu did not converge on iteration step ",
                    iter, "\n")
              }
            }
          }
        }
      }
      # qq = try(optim(c(mu), method = "BFGS", fn = mu_f, 
      #                gr = mu_grad_f, model.coefs = c(new.Lambda, new.B, 
      #                                                new.beta), va.sigma = c(new.sigma), pai = new.pai, 
      #                alpha = new.alpha, control = list(trace = 0, 
      #                                                  fnscale = -1, maxit = maxit, reltol = 1e-06)), 
      #          silent = TRUE)
      # if ("try-error" %in% class(qq)) {
      #   new.mu = mu
      # } else {
      #   if (iter > 1 && m.cur.logfunc > qq$value) {
      #     if (trace) {
      #       cat("Optimization of mu did not improve on iteration step ", 
      #           iter, "\n")
      #     }
      #     new.mu = mu
      #   } else {
      #     if (trace) {
      #       cat("Variational parameters mu updated", 
      #           "\n")
      #     }
      #     new.mu = matrix(qq$par, n, n.factors)
      #     if (qq$convergence != 0) {
      #       if (trace) {
      #         cat("Optimization of mu did not converge on iteration step ", 
      #             iter, "\n")
      #       }
      #     }
      #   }
      # }
      # 
      #if (!is.null(init)) new.mu=init$mu
      q_m = list(value = mu_f(c(new.Lambda, new.B, new.beta), 
                              va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, 
                              alpha = new.alpha))
      m.new.logfunc = q_m$value
      m.cur.logfunc = m.new.logfunc
      q_s = list(value = va.sigma_f(c(new.Lambda, new.B, 
                                      new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), 
                                    pai = new.pai, alpha = new.alpha))
      s.new.logfunc = q_s$value
      s.cur.logfunc = s.new.logfunc
      delta.alpha = abs(new.alpha - alpha)
      sigma = new.sigma
      mu = new.mu
      alpha = new.alpha
      pai = new.pai
      p.iter = p.iter + 1
    }
    q_b = list(value = model.coefs_f(c(new.Lambda, new.B, 
                                       new.beta), va.mu = c(new.mu), va.sigma = c(new.sigma), 
                                     pai = new.pai, alpha = new.alpha))
    b.new.logfunc = q_b$value
    b.cur.logfunc = b.new.logfunc
    qq = list(value = VLB(c(new.Lambda, new.B, new.beta), 
                          va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, 
                          alpha = new.alpha))
    new.VLB = qq$value
    diff = abs(new.VLB - cur.VLB)
    ratio = abs(new.VLB/cur.VLB)
    if (trace) 
      cat("New VLB: ", new.VLB, "cur VLB: ", cur.VLB, "Ratio of VLB: ", 
          ratio, ". Difference in VLB: ", diff, "\n")
    #if (trace) cat("\nThe current evaluate results: ", evaluate_Ex4(promax(new.Lambda)$loadings, promax(data$Lambda)$loadings))
    cur.VLB = new.VLB
    qq = list(value = log_base(c(new.Lambda, new.B, new.beta), 
                               va.mu = c(new.mu), va.sigma = c(new.sigma), pai = new.pai, 
                               alpha = new.alpha))
    cur.log = qq$value
    EE = log_base2(c(new.Lambda, new.B, new.beta), va.mu = c(new.mu), 
                   va.sigma = c(new.sigma), pai = new.pai, alpha = new.alpha)
    EN2 = 2 * cur.log - 2 * EE
    pai = new.pai
    alpha = new.alpha
    B = new.B
    Lambda = new.Lambda
    beta = new.beta
    mu = new.mu
    sigma = new.sigma
    iter = iter + 1
  }
  if (iter > 99) {
    print("FACUB not converging!")
  }
  ll = matrix(new.B, n, p, byrow = TRUE) + X %*% t(new.beta) + 
    new.mu %*% t(new.Lambda)
  e.mat = ll + 0.5 * (new.sigma) %*% t(new.Lambda^2)
  if (is.null(V)) beta = NULL
  out.list$VLB = cur.VLB
  out.list$EE = EE
  out.list$lob = cur.log
  out.list$EN2 = EN2
  out.list$iter = iter - 1
  out.list$lvs$alpha = alpha
  out.list$lvs$mu = mu
  out.list$lvs$sigma = sigma
  out.list$params$pai = pai
  out.list$params$B = B
  out.list$params$beta = beta
  out.list$params$Lambda = Lambda
  out.list$ll = plogis(ll)
  out.list$e.mat = plogis(e.mat)
  return(out.list)
}



generate_data_csda1 = function(n, type = NULL) {
  p = 6
  k = 2
  Lambda = matrix(c(.9, .8 , .7, .5, rep(0, 5), .6, .7, .8), p, k)
  Phi = matrix(c(1, .5, .5, 1), 2, 2)
  if (type == "standard") {
    Phi = diag(k)
  }
  f = mvrnorm(n, rep(0, k), Phi)
  Theta = diag(p) - diag(diag(Lambda%*%Phi%*%t(Lambda)))
  delta = mvrnorm(n, rep(0, p), Theta)
  xstar = f %*% t(Lambda) + delta
  ind1 = which(xstar < -1.2)
  ind2 = which(xstar >= -1.2 & xstar < 0)
  ind3 = which(xstar >= 0 & xstar < 1.2)
  ind4 = which(xstar >= 1.2)
  x = matrix(0, n, p)
  x[ind1] = 1
  x[ind2] = 2
  x[ind3] = 3
  x[ind4] = 4
  return(list(R = x, f = f, Lambda = Lambda, k = k, p = p))
}

generate_data_csda2 = function(n, type = NULL) {
  p = 15
  k = 3
  Lambda = matrix(c(seq(.4, .8, .1), .3, rep(0, p-1), seq(.8, .4, -.1), rep(0, p-1), seq(.5, .9, .1), .4), p, k)
  Phi = matrix(c(1, .2, .5, .2, 1, .8, .5, .8, 1), k, k)
  if (type == "standard") {
    Phi = diag(k)
  }
  f = mvrnorm(n, rep(0, k), Phi)
  Theta = diag(p) - diag(diag(Lambda%*%Phi%*%t(Lambda)))
  delta = mvrnorm(n, rep(0, p), Theta)
  xstar = f %*% t(Lambda) + delta
  ind1 = which(xstar < -1.2)
  ind2 = which(xstar >= -1.2 & xstar < 0)
  ind3 = which(xstar >= 0 & xstar < 1.2)
  ind4 = which(xstar >= 1.2)
  x = matrix(0, n, p)
  x[ind1] = 1
  x[ind2] = 2
  x[ind3] = 3
  x[ind4] = 4
  return(list(R = x, f = f, Lambda = Lambda, k = k, p = p))
}


generate_data_csda = function(n, p, k) {
  Lambda = matrix(runif(p*d, -2, 2), p, k)
  #Lambda = matrix(rnorm(p*k), p, k)
  #Lambda = matrix(c(.9, .8 , .7, .5, rep(0, 5), .6, .7, .8), p, k)
  # Lambda = matrix(0, p, k)
  # for (kk in 1:k) {
  #   Lambda[(p/k*(kk-1)+1):(p/k*kk), kk] = sample(seq(0.5, 1, 0.1), p/k, replace = TRUE)
  # }
  Phi = matrix(0, k, k)
  diag(Phi) = 1
  f = mvrnorm(n, rep(0, k), Phi)
  # f = matrix(rnorm(n*k), n, k)
  Theta = diag(runif(p, 0, 1))
  #Theta = abs(diag(p) - diag(diag(Lambda%*%Phi%*%t(Lambda))))
  delta = mvrnorm(n, rep(0, p), Theta)
  xstar = f %*% t(Lambda) + delta
  ind1 = which(xstar < -1.2)
  ind2 = which(xstar >= -1.2 & xstar < 0)
  ind3 = which(xstar >= 0 & xstar < 1.2)
  ind4 = which(xstar >= 1.2)
  x = matrix(0, n, p)
  x[ind1] = 1
  x[ind2] = 2
  x[ind3] = 3
  x[ind4] = 4
  return(list(R = x, f = f, Lambda = Lambda, k = k, p = p))
}


hic_new = function(R, V=NULL, m, r) {
  n = nrow(R)
  p = ncol(R)
  aic = bic = hic1 = hic2 = hic3 = numeric(r)
  sn = log(n)/2
  for (w in 1:r) {
    initials = generate_initials(R, V, w, type = "3", 0, 0, 0.1, 0.5)
    re = FAVA_cpp(R, V, m, w, initials)
    aic[w] = -2*re$EE + 2*(w*p + n*w-w^2)     # aic
    bic[w] = -2*re$EE + log(n*p)*(w*p + n*w - w^2)  # bic in jasa(2022) 2*n*w?
    hic1[w] = -2*re$EE + log(n)*(w*p-w^2) + 2*(w*n)  # hic1 in JCGS(2021)
    hic2[w] = -re$EE + sn*(w*p -w^2 +w*n- re$EN2)    # hic2 with sn=log(n)/2 in JASA(2023)
    hic3[w] = -re$EE + (w*p -w^2 + w*n- re$EN2)     # hic3 with sn=1 in JASA(2023)
  }
  
  return(list(aic = aic,bic = bic,hic1=hic1, hic2=hic2, hic3=hic3))
}

hic_KLAnnealing_new = function(R, V=NULL, m, r) {
  n = nrow(R)
  p = ncol(R)
  aic = bic = hic1 = hic2 = hic3 = numeric(r)
  sn = log(n)/2
  for (w in 1:r) {
    initials = generate_initials(R, V, w, type = "3", 0, 0, 0.1, 0.5)
    re = FAVA_KLAnnealing_cpp(R, V, m, w, initials)
    aic[w] = -2*re$EE + 2*(w*p + n*w-w^2)     # aic
    bic[w] = -2*re$EE + log(n*p)*(w*p + n*w - w^2)  # bic in jasa(2022)
    hic1[w] = -2*re$EE + log(n)*(w*p-w^2) + 2*(w*n)  # hic1 in JCGS(2021)
    hic2[w] = -re$EE + sn*(w*p -w^2 +w*n- re$EN2)    # hic2 with sn=log(n)/2 in JASA(2023)
    hic3[w] = -re$EE + (w*p -w^2 + w*n- re$EN2)     # hic3 with sn=1 in JASA(2023)
  }
  
  return(list(aic = aic,bic = bic,hic1=hic1, hic2=hic2, hic3=hic3))
}




library(vegan)
library(MASS)
library(tictoc)
library(parallel)
library(Rcpp)

source("functions.R")
gc()
sourceCpp("utils1.cpp")
sourceCpp("MCUB.cpp")

n = 400
p = 20
d = 2
q = 2
m = 5
max_r = ifelse(d < 5, 5, 8)
pai_param = c(.7, 1)

max.cores = detectCores()
num.cores = ifelse(max.cores >= 30, 30, max.cores)
if (d==6 && n==2400) num.cores=10

tic()
outcome = mclapply(1:100, function(seed) {
  set.seed(seed)
  data = generate_data_FACUB(n, p, q, d, m, pai_param)
  R = data$R
  x = data$x
  
  d_minus_1=d-1
  d_plus_1=d+1
  x_minus_1=matrix(data$x[,1:(q-1)],ncol=1)
  x_plus_1=cbind(data$x, rnorm(p))
  initials1 = generate_initials(R, x, d, 0, 0, 0.1, 0.5)
  initials2 = generate_initials(R, x, d_plus_1, 0, 0, 0.1, 0.5)
  initials3 = generate_initials(R, x, d_minus_1, 0, 0, 0.1, 0.5)
  initials4 = generate_initials(R, x_plus_1, d, 0, 0, 0.1, 0.5)
  initials5 = generate_initials(R, x_minus_1, d, 0, 0, 0.1, 0.5)
  initials6 = generate_initials(R, x_minus_1, d_plus_1, 0, 0, 0.1, 0.5)
  initials7 = generate_initials(R, x_minus_1, d_minus_1, 0, 0, 0.1, 0.5)
  initials8 = generate_initials(R, x_plus_1, d_plus_1, 0, 0, 0.1, 0.5)
  initials9 = generate_initials(R, x_plus_1, d_minus_1, 0, 0, 0.1, 0.5)
  
  res1 = FAVA_cpp(R, x, m, d, initials1)
  res2 = FAVA_cpp(R, x, m, d_plus_1, initials2)
  res3 = FAVA_cpp(R, x, m, d_minus_1, initials3)
  res4 = FAVA_cpp(R, x_plus_1, m, d, initials4)
  res5 = FAVA_cpp(R, x_minus_1, m, d, initials5)
  res6 = FAVA_cpp(R, x_minus_1, m, d_plus_1, initials6)
  res7 = FAVA_cpp(R, x_minus_1, m, d_minus_1, initials7)
  res8 = FAVA_cpp(R, x_plus_1, m, d_plus_1, initials8)
  res9 = FAVA_cpp(R, x_plus_1, m, d_minus_1, initials9)
  
  c(evaluate_FACUB2(data, res1),
    evaluate_FACUB2(data, res2),
    evaluate_FACUB2(data, res3),
    evaluate_FACUB2(data, res4),
    evaluate_FACUB2(data, res5),
    evaluate_FACUB2(data, res6),
    evaluate_FACUB2(data, res7),
    evaluate_FACUB2(data, res8),
    evaluate_FACUB2(data, res9))
  
}, mc.cores = num.cores)

toc()

outcome = matrix(unlist(outcome), nrow = 100, byrow = TRUE)
means = colMeans(outcome)
sds = apply(outcome, 2, sd)

# error2, VCC2, and TCC2 are the recovery of Lambda only 
COLnames = c("error1", "VCC1", "TCC1", "error2", "VCC2", "TCC2", "MSE_pai", "MSE_B", "MSE_beta", "MSE_emat", "MSE_ll")
ROWnames = c("d2,q2", "d3,q2", "d1,q2", "d2,q3", "d2,q1", "d3,q1", "d1,q1", "d3,q3", "d1,q3")
means_out = matrix(means, 9, byrow = TRUE, dimnames = list(ROWnames, COLnames))
sds_out = matrix(sds, 9, byrow = TRUE, dimnames = list(ROWnames, COLnames))

df1 = data.frame(means_out)
colnames(df1) = colnames(means_out)
df2 = data.frame(sds_out)
colnames(df2) = colnames(sds_out)

cat("################## Ex2 ###################")
cat("\npai ~ U[", pai_param[1], ",", pai_param[2], "]", ", ")
cat("\nn =", n, ", p =", p, ", q =", q, ", d =", d, ", m =", m)
cat("\nThe means of params: ", "\n")
knitr::kable(df1, digits = 3)
cat("\n", "The sds of params", "\n")
knitr::kable(df2, digits = 3)


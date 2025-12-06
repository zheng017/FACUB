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
p = 10
d = 2
q = 0
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
  
  initials1 = generate_initials(R, x, d, B0 = 0, beta0 = 0, sigma0 = 0.1, pai0 = 0.5)
  
  res1 = FAVA_cpp(R, x, m, d, initials1)
  res2 = MCUB_cpp(R, x, m, b0 = rep(1, p))
  res3 = MCUB_cpp(R, cbind(x, data$f), m, b0 = rep(1, p))
  
  hic1 = hic(R, x, m, max_r)
  
  c(evaluate_FACUB(data, res1),
    evaluate_MCUB(data, res2),
    evaluate_MCUBo(data, res3),
    sapply(hic1, which.min))
}, mc.cores = num.cores)

toc()

outcome = matrix(unlist(outcome), nrow = 100, byrow = TRUE)
means = colMeans(outcome)
sds = apply(outcome, 2, sd)

COLnames = c("error", "VCC", "TCC", "MSE_pai", "MSE_B", "MSE_beta", "MSE_emat", "MSE_ll")
means_out = matrix(means[1:24], 3, byrow = TRUE, dimnames = list(c("FACUB", "MCUB", "MCUB-o"), COLnames))
sds_out = matrix(sds[1:24], 3, byrow = TRUE, dimnames = list(c("FACUB", "MCUB", "MCUB-o"), COLnames))

df1 = data.frame(means_out)
colnames(df1) = colnames(means_out)
df2 = data.frame(sds_out)
colnames(df2) = colnames(sds_out)

cat("################## Ex1 ###################")
cat("\npai ~ U[", pai_param[1], ",", pai_param[2], "]", ", ")
cat("\nn =", n, ", p =", p, ", q =", q, ", d =", d, ", m =", m)
cat("\nThe means of params: ", "\n")
knitr::kable(df1, digits = 3)
cat("\n", "The sds of params", "\n")
knitr::kable(df2, digits = 3)

hic_name = c("aic", "hic1", "hic2")
for (i in 1:3) {
  cat("\n", hic_name[i], ":")
  print(table(outcome[, 24+i]))
}

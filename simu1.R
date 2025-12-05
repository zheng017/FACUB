n=200
p=10
d=2


source("/home/zhengzhi/FACUB/CAUB_func.R")
library(Rcpp)
gc()
sourceCpp("/home/zhengzhi/FACUB/utils1.cpp")
sourceCpp("/home/zhengzhi/FACUB/MCUB.cpp")
# sourceCpp("Rcodes/JML/utils.cpp")
# source("Rcodes/JML/functions.R")
# Ex1: d=1, q=0
library(vegan)
library(MASS)
library(parallel)
library(tictoc)
q = 0
m = 5
#b0 = 1
max_r = ifelse(d < 5, 5, 8)
pai_param = c(.7, 1)
type1 = "uniform"
type2 = "identity"

max.cores = detectCores()
num.cores = ifelse(max.cores >= 30, 30, max.cores)
if (d==6 && n==2400) num.cores=10

tic()
outcome = mclapply(1:100, function(seed) {
  set.seed(seed)
  #data = generate_data(n, p, q, d, m)
  data = generate_data_FACUB(n, p, q, d, m, b0, pai_param, type1, type2)
  R = data$R
  x = data$x
  
  initials1 = generate_initials(R, x, d, type = "3", 0, 0, 0.1, 0.5)
  
  res1 = FAVA_cpp(R, x, m, d, initials1)
  res2 = MCUB_cpp(R, x, m, b0 = rep(1, p))
  res3 = MCUB_cpp(R, cbind(x, data$f), m, b0 = rep(1, p))
  
  hic1 = hic_new(R, x, m, max_r)
  
  c(evaluate_FACUB(data, res1),
    evaluate_MCUB(data, res2),
    evaluate_MCUBo(data, res3),
    sapply(hic1, which.min))
}, mc.cores = num.cores)

toc()

#save.image(paste0(format(Sys.Date(), "%m_%d"), "_", "pai0.7-1_", type1, "_", type2, "_n", n, "_p", p, "_q", q, "_d", d, "_m", m, ".Rdata"))

#outcome = matrix(unlist(outcome), ncol = 17, byrow = TRUE)
outcome = matrix(unlist(outcome), nrow = 100, byrow = TRUE)
means = colMeans(outcome)
sds = apply(outcome, 2, sd)

COLnames = c("error1", "error2", "VCC", "TCC", "MSE_pai", "MSE_B", "MSE_beta", "MSE_emat", "MSE_ll")
means_out = matrix(means[1:27], 3, byrow = TRUE, dimnames = list(c("FACAUB", "MCUB", "MCUB-o"), COLnames))
sds_out = matrix(sds[1:27], 3, byrow = TRUE, dimnames = list(c("FACAUB", "MCUB", "MCUB-o"), COLnames))

df1 = data.frame(means_out)
colnames(df1) = colnames(means_out)
df2 = data.frame(sds_out)
colnames(df2) = colnames(sds_out)

cat("################## Ex1 ###################")
cat("\ntype1 =", type1, ", type2 =", type2, ", pai ~ U[", pai_param2[1], ",", pai_param2[2], "]", ", b0 =", b0, ", ")
cat("\nn =", n, ", p =", p, ", q =", q, ", d =", d, ", m =", m)
cat("\nThe means of params: ", "\n")
knitr::kable(df1, digits = 3)
cat("\n", "The sds of params", "\n")
knitr::kable(df2, digits = 3)

hic_name = c("aic","bic","hic1", "hic2", "hic3")
for (i in 1:5) {
  cat("\n", hic_name[i], ":")
  print(table(outcome[, 27+i]))
}


save.image(paste0(format(Sys.Date(), "%m_%d"), "_", "JIC_", type1, "_", type2, "_n", n, "_p", p, "_q", q, "_d", d, "_m", m, ".Rdata"))


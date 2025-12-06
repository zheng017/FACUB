library(vegan)
library(MASS)
library(lavaan)
library(tictoc)
library(parallel)
library(Rcpp)

source("functions.R")
gc()
sourceCpp("utils1.cpp")
sourceCpp("MCUB.cpp")

n = 200
p = 30
d = 3
m = 4

max_r = ifelse(d < 5, 5, 8)

max.cores = detectCores()
num.cores = ifelse(max.cores >= 30, 30, max.cores)
if (d==6 && n==2400) num.cores=10

tic()
outcome = mclapply(1:100, function(seed) {
  set.seed(seed)
  
  data = generate_data_URV(n, p, d)
  R = data$R
  x = data$x
  
  initials1 = generate_initials(R, x, d, B0 = 0, beta0 = 0, sigma0 = 0.1, pai0 = 0.5)
  
  res1 = FAVA_cpp(R, x, m, d, initials1)
  res2 = efa(data.frame(R), nfactors = d, ordered = TRUE)
  
  hic1 = hic(R, x, m, max_r)
  
  c(evaluate_recovery(data$Lambda, res1$params$Lambda),
    evaluate_recovery(data$Lambda, res2$loadings),
    sapply(hic1, which.min))
}, mc.cores = num.cores)

toc()

outcome = matrix(unlist(outcome), nrow = 100, byrow = TRUE)
means = colMeans(outcome)
sds = apply(outcome, 2, sd)

means_out = matrix(means[1:6], 2, byrow = TRUE, dimnames = list(c("FACAUB", "lavaan"), c("error", "VCC", "TCC")))
sds_out = matrix(sds[1:6], 2, byrow = TRUE, dimnames = list(c("FACAUB", "lavaan"), c("error", "VCC", "TCC")))

df1 = data.frame(means_out)
colnames(df1) = colnames(means_out)
df2 = data.frame(sds_out)
colnames(df2) = colnames(sds_out)

cat("################## Ex3 ###################")
cat("\nn =", n, ", p =", p, ", d =", d, ", m =", m)
cat("\nThe means of params: ", "\n")
knitr::kable(df1, digits = 3)
cat("\n", "The sds of params", "\n")
knitr::kable(df2, digits = 3)


hic_name = c("aic", "hic1", "hic2")
for (i in 1:3) {
  cat("\n", hic_name[i], ":")
  print(table(outcome[, 6+i]))
}

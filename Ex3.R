library(vegan)
library(MASS)
library(lavaan)
library(tictoc)
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


data = generate_data_URV(n, p, d)
R = data$R
x = data$x

initials1 = generate_initials(R, x, d, B0 = 0, beta0 = 0, sigma0 = 0.1, pai0 = 0.5)

res1 = FAVA_cpp(R, x, m, d, initials1)
res2 = efa(data.frame(R), nfactors = d, ordered = TRUE)

hic1 = hic(R, x, m, max_r)

evaluate_recovery(data$Lambda, res1$params$Lambda)
evaluate_recovery(data$Lambda, res2$loadings)
sapply(hic1, which.min)

library(vegan)
library(MASS)
library(tictoc)
library(Rcpp)

source("functions.R")
gc()
sourceCpp("utils1.cpp")
sourceCpp("MCUB.cpp")

n = 400
p = 20
d = 2
q = 0
m = 5
max_r = ifelse(d < 5, 5, 8)
pai_param = c(.7, 1)

data = generate_data_FACUB(n, p, q, d, m, pai_param)
R = data$R
x = data$x

initials1 = generate_initials(R, x, d, B0 = 0, beta0 = 0, sigma0 = 0.1, pai0 = 0.5)

res1 = FAVA_cpp(R, x, m, d, initials1)
res2 = MCUB_cpp(R, x, m, b0 = rep(1, p))
res3 = MCUB_cpp(R, cbind(x, data$f), m, b0 = rep(1, p))

hic1 = hic(R, x, m, max_r)

evaluate_FACUB(data, res1)
evaluate_MCUB(data, res2)
evaluate_MCUBo(data, res3)
sapply(hic1, which.min)

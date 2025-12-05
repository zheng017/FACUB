n=200
p=10
d=2

source("CAUB_func.R")
library(Rcpp)
gc()
sourceCpp("/home/zhengzhi/FACUB/utils1.cpp")
sourceCpp("/home/zhengzhi/FACUB/MCUB.cpp")

library(vegan)
library(MASS)
library(tictoc)
q = 0
m = 5
#b0 = 1
max_r = ifelse(d < 5, 5, 8)
pai_param = c(.7, 1)
type1 = "uniform"
type2 = "identity"

#set.seed(seed)
#data = generate_data(n, p, q, d, m)
data = generate_data_FACUB(n, p, q, d, m, b0, pai_param, type1, type2)
R = data$R
x = data$x

initials1 = generate_initials(R, x, d, type = "3", 0, 0, 0.1, 0.5)

res1 = FAVA_cpp(R, x, m, d, initials1)
res2 = MCUB_cpp(R, x, m, b0 = rep(1, p))
res3 = MCUB_cpp(R, cbind(x, data$f), m, b0 = rep(1, p))

hic1 = hic_new(R, x, m, max_r)

evaluate_FACUB(data, res1
evaluate_MCUB(data, res2)
evaluate_MCUBo(data, res3)
sapply(hic1, which.min))


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
q = 2
m = 5
max_r = ifelse(d < 5, 5, 8)
pai_param = c(.7, 1)

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

evaluate_FACUB2(data, res1)
evaluate_FACUB2(data, res2)
evaluate_FACUB2(data, res3)
evaluate_FACUB2(data, res4)
evaluate_FACUB2(data, res5)
evaluate_FACUB2(data, res6)
evaluate_FACUB2(data, res7)
evaluate_FACUB2(data, res8)
evaluate_FACUB2(data, res9)
  

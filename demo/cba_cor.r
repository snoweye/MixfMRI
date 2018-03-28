library(MixfMRI, quietly = TRUE)
dim <- c(4, 5, 4, 10)
set.seed(123)
da.ts <- array(rnorm(prod(dim)), dim = dim)
id.class <- cba.cor(da.ts)
table(id.class)

fun.tanh <- function(a, B){
  d <- 1 / apply(B, 2, function(b){ dist(rbind(as.vector(a), b)) })
  tanh(d)
}
id.class.tanh <- cba.cor(da.ts, fun.sim = fun.tanh)
table(id.class.tanh)

fun.logit <- function(a, B){
  d <- dist(t(cbind(a, B)))[1:ncol(B)]
  (1 / (1 + exp(-d))) * 2 - 1
}
id.class.logit <- cba.cor(da.ts, fun.sim = fun.logit)
table(id.class.logit)


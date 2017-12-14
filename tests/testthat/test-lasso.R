library(microbenchmark)
library(glmnet)
set.seed(749125)
context("Unit test for the generic LASSO estimation procedure.")

# generate simple test data
n = 1000
p = 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean=0, sd=0.1)

testn <- 1e4
testx <- matrix(rnorm(testn * p), testn, p)
testy <- sin(testx[, 1]) * sin(testx[, 2]) + rnorm(testn, mean=0, sd=0.1)

# fit design matrix for HAL
basis_list <- hal9001:::enumerate_basis(x)
x_basis <- hal9001:::make_design_matrix(x, basis_list)
time_design_matrix <- proc.time()

# catalog and eliminate duplicates
copy_map <- hal9001:::make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]


#################################################
# test generating the sequence of lambdas



xcenter <- hal9001:::get_pnz(x_basis)
xscale <- hal9001:::get_xscale(x_basis)
ybar <- mean(y)
y_centered <- y - ybar
# xscale = rep(1, length(xscale))
lambda_max <- hal9001:::find_lambda_max(x_basis, y_centered, xscale, xcenter)

# verify that lambda max zeros out coefs
beta <- rep(0, ncol(x_basis))
lassi_1step <- hal9001:::lassi_fit_cd(X = x_basis, resids = y_centered,
                                      beta = beta, lambda = 0,
                                      nsteps = 1, xscale = xscale, xcenter=xcenter,
                                      intercept=ybar, active_set = FALSE)
beta
lassi_1step <- hal9001:::lassi_fit_cd(X = x_basis, resids = y_centered,
                                      beta = beta, lambda = lambda_max,
                                      nsteps = 1, xscale = xscale, xcenter=xcenter,
                                      intercept=ybar, active_set = FALSE)
test_that("lambda_max results in zero beta vector", expect_true(all(beta == 0)))

# verify that a slightly smaller lambda does not
delta <- 1 - 1e-3
lambda_delta <- lambda_max * delta
beta0 <- rep(0, ncol(x_basis))
lassi_smallstep <- hal9001:::lassi_fit_cd(X = x_basis, resids = y_centered,
                                          beta = beta, lambda = lambda_delta,
                                          nsteps = 1, xscale = xscale, xcenter=xcenter,
                                          intercept=ybar, active_set = FALSE)
test_that("a slightly smaller lambda results in nonzero beta vector",
          expect_true(!all(beta == 0)))

# generate sequence of lambdas
nlambda <- 100
lambda_min_ratio <- 0.01
lambdas <- hal9001:::lambda_seq(lambda_max = lambda_max,
                                lambda_min_ratio = lambda_min_ratio,
                                nlambda = nlambda)
test_that("lambda_seq generates a sequence of lambdas",{
  expect_length(lambdas, nlambda)
  expect_equal(max(lambdas), lambda_max)
  expect_equal(min(lambdas), lambda_max * lambda_min_ratio)
})

#################################################
# test a single coordinate descent update
ybar=mean(y)
resid <- y[seq_along(y)] - ybar
resid2 <- y[seq_along(y)] - ybar
beta <- rep(0, ncol(x_basis))
pre_mse <- mean(resid^2)

n <- length(y)
i <- 123 #which beta to update (1 - indexed)

xscale <- hal9001:::get_xscale(x_basis)
xcenter <- hal9001:::get_pnz(x_basis)

#explicitly scale for verification methods
xvar <- x_basis[, i] / xscale[i]
xv2 <- mean((xvar-mean(xvar))^2)
test_that("xscale correctly scales x_basis", expect_equal(xv2, 1))
cp=X_t_resid(x_basis, resid, i-1, xscale[i], xcenter[i])/n

get_new_beta(x_basis, resid, i-1, 0, xscale[i], xcenter[i])
ls_beta <- coef(lm(resid ~ xvar - 1))
cp <- as.vector(resid%*%x_basis[,i]/xscale[i])
cd_beta <-  cp / n;
ls_beta/cd_beta
sum(x_basis[,i])

xscale[i]^2/(1-xcenter[i])
coord_update <- hal9001:::update_coord(x_basis, resid, beta, 0, i - 1, xscale, xcenter)
beta_new <- beta[i]
test_that("coordinate descent works as it does in R", expect_equal(cd_beta,
                                                                   beta_new))

post_mse <-  mean(resid^2)
verify_resid = y - (ybar + beta[i] * x_basis[,i]/xscale[i])
all.equal(resid, verify_resid)

verify_postmse  <-  mean(verify_resid^2)
test_that("the mse of the updated residuals is as expected",
          expect_equal(post_mse, verify_postmse))


microbenchmark({
  glmnet::glmnet(x = x_basis, y = y, intercept=TRUE, nlambda = 100,
  lambda.min.ratio = 0.01, family = "gaussian", alpha = 1)
})
microbenchmark({
  lassi(x_basis, y, nlambda=100, lambda_min_ratio = 0.01)
}, times = 1)

################################################################################
# PREDICTION
################################################################################

# format test data set
new_data <- as.matrix(testx)
pred_x_basis <- hal9001:::make_design_matrix(new_data, basis_list)

pred_x_basis <- hal9001:::apply_copy_map(pred_x_basis, copy_map)
system.time({
# lassi prediction and mses
lassi_fit <- hal9001:::lassi(x_basis, y)
})

system.time({
beta <- lassi_fit$beta_mat
xscale[xscale==0]=min(xscale[xscale!=0])
beta <- beta / xscale
# intercept <- lassi_fit$intercept - crossprod(xcenter, beta)
lassi_fit$beta_mat <- beta
# lassi_fit$intercept <- intercept
})
pred_mat <- predict(lassi_fit, pred_x_basis)
mses <- apply(pred_mat, 2, function(preds) {mean((preds - testy)^2)})


# glmnet predictions and mses
g <- glmnet::glmnet(x = x_basis, y = y, intercept = TRUE,
                    nlambda = 100, lambda.min.ratio = 0.01, family = "gaussian",
                    alpha = 1, standardize.response = FALSE, standardize = TRUE)
glmnet_beta_mat <- coef(g)
test_that("fits use same lambda vector", expect_equal(g$lambda, lassi_fit$lambdas))

pred_mat <- predict(g, pred_x_basis)
gmses <- apply(pred_mat, 2, function(preds) {mean((preds - testy)^2)})

test_that("lassi isn't doing much worse in terms of MSE",
          expect_lt(max((mses-gmses)/gmses), 1e-8))

li = lassi_fit$intercept
gi = g$a0
lb = lassi_fit$beta_mat
gb = g$beta

plot(li, gi)
library(biglasso)
microbenchmark({
xbm <- as.big.matrix(as.matrix(x_basis))
}, times=1)
microbenchmark({
bl <- biglasso(xbm, y, lambda=g$lambda, screen="None")
}, times=1)

blcoef=coef(bl)[-1,]
gcoef=coef(g)[-1,]
plot(blcoef[,12], gcoef[,12])
plot(apply(abs(lcoef-blcoef),2,max))
lcoef=lassi_fit$beta_mat
pxbm <- as.big.matrix(as.matrix(pred_x_basis))

blpred=predict(bl, pxbm)
blmses <- apply(blpred, 2, function(preds) {mean((preds - testy)^2)})

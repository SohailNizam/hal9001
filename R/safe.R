# lasso_safe=function(){
#   dot=(resid%*%x)/(xscale*n)
#   ynorm=norm(resid,type="2")/n
#   rho=as.vector((ynorm+abs(dot))/(ynorm+lambda_max))
#   thresh=rho*lambda_max
#   betas=lasso_init$beta_mat
#   alambda=lambdas[1]
#   lambda_min=min(lambdas)
#   sum(alambda>thresh) #num zero betas according to safe
#   nnz_betas=apply(betas!=0,2,sum)
# }
# 
# #so what we want to do, is go through and find the nz coefs at each step, and then only update those
# #what we're doing now is to update nz coefs, then update all coefs (inefficient)
# #basically, given the previous solution, we can examine the gradient and compare it to the new lambda value
# #if it's below some threshold, we know there's no update
# 
# #this isn't too different from soft thresholding in the first place
# #so safe specifies a lambda threshold above which betas are _guaranteed_ to be zero
# #it does not make any claims about betas below that threshold (ie they could still be zero)
# #so if we wanted to solve a lambda problem with dimension, say 107, w
# sort(thresh)
# lambda5=0.12587
# lambda6=0.12015
# betas=c(357,571,349)
# #so rho[349]*lambda_max should be between lambda5 and lambda6
# #so rho[349] should be around 
# ideal_rho=lambdas[6]/lambda_max
# dot349=dot[349]
# #so we want to solve ideal_rho=(a+dot[349])/(a+lambda_max)
# # ideal_rho*lambda_max-dot[349]=a*(1-ideal_rho)
# # a=(ideal_rho*lambda_max-dot[349])/(1-ideal_rho)
# a=(ideal_rho*lambda_max-dot[349])/(1-ideal_rho)
# test_rho=(a+dot[349])/(a+lambda_max)
# beta_step=apply(betas!=0,1,function(x)min(which(x)))
# rev(order(dot))[1:100]
# order(beta_step)[1:100]
# a=rank(beta_step)
# b=rank(max(dot)-abs(dot))
# plot(a,b)
# #guaranteed to be absent (maybe finds a subset of all zeros?)
# thresh=lambdas-(ynorm/n)*(lambda_max-lambdas)/lambda_max
# lambdas
# lambda_max

# reisd
# resid
# microbenchmark({
# grad=resid%*%x_basis
# })

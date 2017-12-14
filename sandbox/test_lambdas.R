# n = 100
# x = cbind(rnorm(n, 1, 0.5),rnorm(n, 1, 0.5))
# y = 2*x[,1]+rnorm(n, 0, 0.1)
x=x_basis
g=glmnet(x=x, y=y)
glmnet_lambdas=g$lambda
var2=function(x){mean((x-mean(x))^2)}
xscale=sqrt(apply(x,2,var2))
xscale2=get_xscale(x)
sum(abs(xscale-xscale2)<1e-4)==ncol(x)-1
scaled=sweep(as.matrix(x_basis),2, xscale, `/`)
scale_check=apply(scaled,2, var2)
max(abs(scale_check-1),na.rm=T)<1e-10

scaled2=sweep(as.matrix(x_basis),2, xscale2, `/`)
scale_check2=apply(scaled2,2, var2)
max(abs(scale_check2[-1000]),na.rm=T)<1e-2

table(scale_check2, exclude=c())
dotprod=as.vector((y-mean(y))%*%x/xscale)
lambda_guess=max(abs(dotprod[!is.nan(dotprod)]))/n
print(lambda_guess[!is.nan(lambda_guess)])
print(glmnet_lambdas[1])

plot(xscale2,xscale)

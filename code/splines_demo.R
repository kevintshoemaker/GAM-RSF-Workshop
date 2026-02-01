# Splines demo for GAM workshop


rm(list=ls())

## (examples modified from Simon Wood's GAM book)

library(gamair)
data(engine)
?engine   # learn more about dataset
engine = engine[order(engine$size),]  # sort by cylinder capacity in liters
size=engine$size; wear=engine$wear   # safe "attach"

# explore orthogonal polynomial basis expansion --------

xseq = seq(min(size),max(size),length=100)
X = cbind(rep(1,100), poly(xseq,5) )    # polynomial basis expansion and intercept

cor(X)   # note that basis functions are uncorrelated

## visualize poly basis -------

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
cols=rainbow(6)
lines(xseq,X[,1],col=cols[1],lwd=2)
lines(xseq,X[,2],col=cols[2],lwd=2)
lines(xseq,X[,3],col=cols[3],lwd=1)
lines(xseq,X[,4],col=cols[4],lwd=1,lty=2)
lines(xseq,X[,5],col=cols[5],lwd=1,lty=3)
lines(xseq,X[,6],col=cols[6],lwd=1,lty=4)

## weighted polynomial basis expansion -----

w = c(3,-1.5,0.6,-0.5,-0.2,-0.1)

lp = X %*% w   # the weird symbols are matrix multiplication

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
lines(xseq,X[,1]*w[1],col=cols[1],lwd=2)
lines(xseq,X[,2]*w[2],col=cols[2],lwd=2)
lines(xseq,X[,3]*w[3],col=cols[3],lwd=1)
lines(xseq,X[,4]*w[4],col=cols[4],lwd=1,lty=2)
lines(xseq,X[,5]*w[5],col=cols[5],lwd=1,lty=3)
lines(xseq,X[,6]*w[6],col=cols[6],lwd=1,lty=4)

lines(xseq,lp[,1],lwd=2,col="darkgreen")  # 

## use lm to find optimal weights ------

Xo = cbind(rep(1,length(size)), poly(size,5) )    # polynomial basis expansion and intercept
m = lm(wear~0+Xo)
w_fit = coef(m)

lp = X %*% w_fit   # fitted linear predictor from lm()

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
lines(xseq,lp[,1],lwd=2,col="black")  # 

coef(m)   # optimal weights determined by 'lm'

## explore other possible functional forms --------

     # random weights- generate possible functional forms
ws = replicate(10,(X%*%runif(ncol(X),-2,2))[,1] )

plot(1,1,xlim=c(1.3,3.1),ylim=c(-5,6),pty="n")
sapply(1:ncol(ws), function(t) lines(xseq,ws[,t],col=sample(cols,1))  )


# explore 'tent' basis ----------
   ## also known as 'piecewise linear regression'

j=3

# tent function: define the basis associated with each "knot" j
tf <- function(x,xj,j){  # x is some arbitrary covariate vector, xj is vector of knot locs, j is the jth knot
  dj=xj*0; dj[j]=1  # y vector is all zero except for the focal knot
  approx(xj,dj,x)$y   # jth tent basis function- fades away as get further from focal knot
}

# aside- how 'approx' function generates a tent basis
covar = runif(100,2,4)
nknots = 10
thisknot = 3  # generate third knot
knots=seq(min(covar),max(covar),length=nknots)   # define knot locations that span the covariate values
y=numeric(nknots); y[thisknot]=1   # focal knot is 0, all others are 1
bf=approx(knots,y,covar)$y  # use approx function to generate a 'tent' that represents distance to focal knot
plot(covar,bf,type="p")

# function for generating model matrix representing a "tent function" basis
tf.X <- function(x,xj){
  nk = length(xj); n=length(x)
  X = matrix(NA,n,nk)
  for(j in 1:nk) X[,j] <- tf(x,xj,j)
  X   # return the basis
}

sj=seq(min(size),max(size),length=6)  # define knots

X = tf.X(xseq,sj)  # design matrix

## visualize tent basis -------

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
cols=rainbow(6)
lines(xseq,X[,1],col=cols[1],lwd=2)
lines(xseq,X[,2],col=cols[2],lwd=2)
lines(xseq,X[,3],col=cols[3],lwd=1)
lines(xseq,X[,4],col=cols[4],lwd=1,lty=2)
lines(xseq,X[,5],col=cols[5],lwd=1,lty=3)
lines(xseq,X[,6],col=cols[6],lwd=1,lty=4)

## visualize range of functional forms -----------

ws = replicate(10,(X%*%runif(ncol(X),-3,3))[,1] )

plot(1,1,xlim=c(1.3,3.1),ylim=c(-5,6),pty="n")
sapply(1:ncol(ws), function(t) lines(xseq,ws[,t],col=sample(cols,1))  )


## find optimal weights using lm -------

Xo = tf.X(size,sj)
b = lm(wear~0+Xo)  # piecewise linear regression

plot(size,wear)
lines(xseq,X %*% coef(b),lwd=2)


## penalized smooth fitting --------

y=wear;x=size;sp=2; xj=sj  # sp is complexity penalty
prs.fit <- function(y,x,xj,sp){
  X=tf.X(x,xj)  # model matrix (tent basis)
  D = diff(diag(length(xj)),differences=2) # columnwise second difference: difference of difference! (analog to second derivative)
  X=rbind(X,sqrt(sp)*D)    # add the penalty
  y=c(y,rep(0,nrow(D)))
  lm(y~0+X)
}

sj=seq(min(size),max(size),length=20)  # note: now we have 20 basis functions
b=prs.fit(y=wear,x=size,xj=sj,sp=5)   ## try different smoothing penalties
plot(size,wear)
Xp = tf.X(xseq,sj)
lines(xseq,Xp%*% coef(b))

# cubic regression splines -------

library(splines)

X = splines::bs(xseq,df=6,intercept=T)   # design matrix for cubic regression splines

## visualize cr basis -------

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
cols=rainbow(6)
lines(xseq,X[,1],col=cols[1],lwd=2)
lines(xseq,X[,2],col=cols[2],lwd=2)
lines(xseq,X[,3],col=cols[3],lwd=1)
lines(xseq,X[,4],col=cols[4],lwd=1,lty=2)
lines(xseq,X[,5],col=cols[5],lwd=1,lty=3)
lines(xseq,X[,6],col=cols[6],lwd=1,lty=4)


## visualize range of functional forms -----------

ws = replicate(10,(X%*%runif(ncol(X),-3,3))[,1] )

plot(1,1,xlim=c(1.3,3.1),ylim=c(-5,6),pty="n")
sapply(1:ncol(ws), function(t) lines(xseq,ws[,t],col=sample(cols,1))  )


## find optimal weights using lm -------

Xo = bs(size,df=6,intercept = T)
b = lm(wear~0+Xo)  # piecewise linear regression

plot(size,wear)
lines(xseq,X %*% coef(b),lwd=2)


## penalized smooth fitting --------

library(mgcv)
mod = gam(wear~s(size,bs="cr",k=8,sp=1),method="REML")
b=coef(mod)   ## try different smoothing penalties
plot(size,wear)
Xp= predict(mod,newdata = data.frame(size=xseq), type = "lpmatrix")
lines(xseq,Xp%*%b, lwd=2)

## try GAM defaults
mod = gam(wear~s(size,bs="cr",k=8),method="REML")
b=coef(mod)   ## try different smoothing penalties
plot(size,wear)
Xp= predict(mod,newdata = data.frame(size=xseq), type = "lpmatrix")
lines(xseq,Xp%*%b, lwd=2)


## END SCRIPT 





















# Portuguese Lark example (presence absence data)
   # example from Simon Wood's GAM book

rm(list=ls())

# load packages -------

library('gamair')
library(ggplot2)
library(mgcv)
library(tidyverse)
library(gratia)


# load data -------------

data(bird)

# logistic regression with GAM -----------

## process data --------

bird <- transform(bird,
                  crestlark = factor(crestlark),
                  linnet = factor(linnet),
                  e = x / 1000,
                  n = y / 1000)
head(bird)

## visualize data --------

ggplot(bird, aes(x = e, y = n, colour = crestlark)) + 
  geom_point(size = 0.5) + coord_fixed() + 
  scale_colour_discrete(na.value = '#bbbbbb33') + 
  labs(x = NULL, y = NULL)


## run model ------

   ## smoother is used to account for spatial autocorrelation
       ## soak up variance due to unmeasured spatial covariates
       ## similar to kriging- interpolate value at all locations based on spatial correlation

bird = subset(bird,!is.na(crestlark))
crest <- gam(crestlark ~ s(e, n, k = 100),
             data = bird,
             family = binomial,
             method = 'REML')

summary(crest)

plot(crest)
draw(crest,rug = F)


## model checking -----

library(DHARMa)

res=simulateResiduals(crest)
plot(res)
plotResiduals(res)
DHARMa::testResiduals(res) # battery of tests

DHARMa::plotResiduals(res)
DHARMa::plotQQunif(res)

DHARMa::testSpatialAutocorrelation(res,x=bird$e,y=bird$n)

DHARMa::testUniformity(res)
DHARMa::testDispersion(res)
DHARMa::testOutliers(res)
# DHARMa::testOutliers(res, type="bootstrap")
DHARMa::testQuantiles(res)

DHARMa::testZeroInflation(res)

# par(mfrow = c(1,2))
plotResiduals(res, bird$e)
plotResiduals(res, bird$n)

plotResiduals(res, bird$e, quantreg = T)
plotResiduals(res, bird$n, quantreg = T)

# NOTE: generally, binomial 0/1 response is difficult to assess with diagnostic tests
  # can be more informative to aggregate over larger areas to obtain a count variable


# Binomial (versus Bernoulli/logistic) regression ------

## aggregate to 10 km scale ------

## convert back to numeric
    # note: transform and mutate are doing basically the same thing
bird <- transform(bird,
                  crestlark = as.numeric(as.character(crestlark)),
                  linnet = as.numeric(as.character(linnet)))

## some variables to help aggregation
bird <- transform(bird, tet.n = rep(1, nrow(bird)),
                  N = rep(1, nrow(bird)), stringsAsFactors = FALSE)

## set to NA if not surveyed
bird$N[is.na(as.vector(bird$crestlark))] <- NA

## aggregate by "quadricula"

bird2 <- bird %>% 
  group_by(QUADRICULA) %>% 
  summarize(
    lark=sum(crestlark,na.rm=T),
    nsurv = sum(N,na.rm=T),
    east = mean(e),
    north = mean(n)
  )

table(bird2$nsurv)

bird2 <- bird2[bird2$nsurv>0,]  # remove zero effort grid cells



## fit aggregated binomial GAM  -------

crest2 <- gam(cbind(lark, nsurv - lark) ~ s(east, north, k = 100),
              data = bird2, family = binomial, method = 'REML')

draw(crest2,rug=F)

summary(crest2)

k.check(crest2)

## model checking -----

ggplot(data.frame(Fitted = fitted(crest2),
                  Resid = resid(crest2)),
       aes(Fitted, Resid)) + geom_point() 


library(DHARMa)

res=simulateResiduals(crest2)
testResiduals(res)

testZeroInflation(res)  # evidence for zero inflation?
testSpatialAutocorrelation(res,x=bird2$east,y=bird2$north)

#  try overdispersed binomial

crest3 <- gam(cbind(lark, nsurv - lark) ~ s(east, north, k = 100),
              data = bird2, family = quasibinomial(), method = 'REML')

draw(crest3)

res=simulateResiduals(crest3)  # oops- DHARMa doesn't work for quasibinomial


## try zero inflated model --------

## NOTE: zero inflated binomial gave me trouble- much harder than it looks
   # instead, run a zero inflated poisson with offset for effort (num surveys)

crest4 <- gam(lark ~ s(east, north, k = 100) + offset(log(nsurv)),
              data = bird2, family = "ziP", method = 'REML')

draw(crest4)

res=simulateResiduals(crest4)
testResiduals(res)
testZeroInflation(res)      # seems much better


  # try poisson with individual random effects on the poisson
crest5 <- gam(list(lark ~ s(east, north, k = 100) + offset(log(nsurv)),~1),
              family=ziplss(),  data = bird2, method = 'REML')

res=simulateResiduals(crest5)    # DHARMa quantile residuals 

# testResiduals(res)  # doesn't really work with DHARMa test
testOutliers(res) 
testZeroInflation(res)

draw(crest5,rug=F)


# use 'sdmTMB' package, which we will use again later   --------------

library(sdmTMB)

# Scale coordinates (recommended for GP models)
bird2$east_scaled <- bird2$east/100
bird2$north_scaled <- bird2$north/100

mesh <- make_mesh(bird2, xy_cols = c("east_scaled", "north_scaled"), cutoff = .1)

plot(mesh)

# mesh <- make_mesh(bird2, xy_cols = c("east", "north"), cutoff = 10)
# plot(mesh)

fit <- sdmTMB(
  cbind(lark, nsurv - lark) ~ 1,
  data = bird2,
  mesh = mesh,
  family = binomial(),
  spatial = "on"
)

summary(fit)
sanity(fit)

x_range <- range(bird2$east_scaled)
y_range <- range(bird2$north_scaled)

# Create a regular grid
grid_points <- expand.grid(
  east_scaled = seq(x_range[1], x_range[2], length.out = 50), # Adjust length.out for grid density
  north_scaled = seq(y_range[1], y_range[2], length.out = 50)
)

p <- predict(fit, newdata = grid_points)

p$east = p$east_scaled*sd(bird2$east) + mean(bird2$east)
p$north = p$north_scaled*sd(bird2$north) + mean(bird2$north)

# p <- predict(fit)

ggplot(p, aes(east, north, z = est)) + geom_contour_filled() #+
  # scale_fill_viridis_c(trans = "sqrt")

## check fit --------

    # randomized quantile residuals
res = residuals(fit, type = "mle-mvn")
# rq_res <- rq_res[is.finite(rq_res)] # remove Inf if present
qqnorm(res);abline(0, 1)

# mcmc (relaxes normality assumption for random effects)

# install.packages("remotes")     # install sdmTMBextra package if not yet installed
# remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)

    # don't worry if this doesn't run!
# samps <- sdmTMBextra::predict_mle_mcmc(fit, mcmc_iter = 800, mcmc_warmup = 400)
# res <- residuals(fit, type = "mle-mcmc", mcmc_samples = samps)
# qqnorm(res);abline(0, 1)

    ### simulation based quantile residuals and DHARMa!   --------

s <- simulate(fit, nsim = 500, type = "mle-mvn")
 
     # check zero inflation
sum(bird2$lark == 0) / length(bird2$lark)   # 52% zeros (observed)
sum(s == 0)/length(s)    # 46% zeros (simulated)
dharma_residuals(s,fit)  #

   # ability to run other DHARMa tests!
res = dharma_residuals(s, fit, return_DHARMa = TRUE)

testResiduals(res)

testZeroInflation(res)

testSpatialAutocorrelation(res,bird2$east,bird2$north)

plot(res)


## try overdispersed binomial (betabinomial family) ---------

fit <- sdmTMB(
  cbind(lark, nsurv - lark) ~ 1,
  data = bird2,
  mesh = mesh,
  family = betabinomial(),
  spatial = "on"
)

summary(fit)
sanity(fit)

x_range <- range(bird2$east_scaled)
y_range <- range(bird2$north_scaled)

# Create a regular grid
grid_points <- expand.grid(
  east_scaled = seq(x_range[1], x_range[2], length.out = 50), # Adjust length.out for grid density
  north_scaled = seq(y_range[1], y_range[2], length.out = 50)
)

p <- predict(fit, newdata = grid_points)

p$east = p$east_scaled*sd(bird2$east) + mean(bird2$east)
p$north = p$north_scaled*sd(bird2$north) + mean(bird2$north)

# p <- predict(fit)

ggplot(p, aes(east, north, z = est)) + geom_contour_filled() #+

## check fit --------

### simulation based quantile residuals and DHARMa!

s <- simulate(fit, nsim = 500, type = "mle-mvn")

sum(bird2$lark == 0) / length(bird2$lark)   # 52% zeros (observed)
sum(s == 0)/length(s)    # 47% zeros (simulated)
dharma_residuals(s,fit)  #

# ability to run other DHARMa tests!
res = dharma_residuals(s, fit, return_DHARMa = TRUE)

testResiduals(res)

testZeroInflation(res)

testSpatialAutocorrelation(res,bird2$east,bird2$north)

### looks better...


# END SCRIPT



### BONUS

    # penalized complexity prior on range and sigma
?pc_matern
pc <- pc_matern(range_gt = 0.01, sigma_lt = 10)

fit <- sdmTMB(
  cbind(lark, nsurv - lark) ~ 1,
  data = bird2,
  mesh = mesh,
  priors = sdmTMBpriors(matern_s = pc),   # note the prior!
  family = binomial(),
  spatial = "on"
)

p <- predict(fit, newdata = grid_points)

p$east = p$east_scaled*sd(bird2$east) + mean(bird2$east)
p$north = p$north_scaled*sd(bird2$north) + mean(bird2$north)

# p <- predict(fit)

ggplot(p, aes(east, north, z = est)) + geom_contour_filled() #+
# scale_fill_viridis_c(trans = "sqrt")







### Basic example from sdmTMB help files - with additional extensions!

### This script covers all the key features of sdmTMB package

rm(list=ls())

# load packages ------

library(dplyr)
library(ggplot2)
library(sdmTMB)
head(pcod)


# visualize data ------

ggplot(pcod, aes(X, Y, col = density)) +
  geom_point() +
  coord_fixed() +
  scale_colour_viridis_c(trans = "sqrt") +
  labs(colour = "Density", title = "Observed survey tows")

glimpse(pcod)



# make mesh ------

mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)

plot(mesh)


# basic logistic regression -------

m <- sdmTMB(
  data = pcod,
  formula = present ~ depth_scaled + depth_scaled2,
  family = binomial(link = "logit"),
  spatial = "off"
)
m

## compare with glm

m0 <- glm(
  data = pcod,
  formula = present ~ depth_scaled + depth_scaled2,
  family = binomial(link = "logit")
)
summary(m0)

# fit basic model with spatial field --------

m1 <- sdmTMB(
  data = pcod,
  formula = present ~ depth_scaled + depth_scaled2,
  mesh = mesh,
  family = binomial(link = "logit"),
  spatial = "on"
)
m1

## try density as response --------
fit <- sdmTMB(
  density ~ s(depth),   # GAM-style smoother on depth
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on"
)

fit

## explore results --------
tidy(fit, conf.int = TRUE)


## check results -------
sanity(fit)

## plot partial effects -----

ggeffects::ggpredict(fit, "depth [50:400, by=2]") |> plot()

ggeffects::ggpredict(fit, "depth [0:600, by=2]") |> plot()

visreg::visreg(fit, "depth")   # on link scale...

## predict on new data ------

glimpse(qcs_grid)
p <- predict(fit, newdata = qcs_grid)


## visualize predictions across study area --------

ggplot(p, aes(X, Y, fill = exp(est))) + geom_raster() +
  scale_fill_viridis_c(trans = "sqrt")


plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed()
}

   # nicer map
plot_map(p, exp(est)) +
  scale_fill_viridis_c(
    trans = "sqrt",
    # trim extreme high values to make spatial variation more visible:
    na.value = "yellow", limits = c(0, quantile(exp(p$est), 0.995))
  ) +
  ggtitle("Prediction (fixed effects + all random effects)",
          subtitle = paste("maximum estimated biomass density =", round(max(exp(p$est))))
  )

# try binomial presence absence model ----------

fit <- sdmTMB(
  present ~ s(depth),
  data = pcod, 
  mesh = mesh,
  family = binomial(link = "logit")
)

fit


# try hurdle/delta model ----------

   # takes a little longer...
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  mesh = mesh,
  family = delta_gamma(link1 = "logit", link2 = "log"),
  spatial="on"
)

fit


# try spatiotemporal model ----------

m2 <- sdmTMB(
  data = pcod,
  formula = present ~ depth_scaled + depth_scaled2,
  mesh = mesh,
  family = binomial(link = "logit"),
  spatial = "on",
  time = "year",
  spatiotemporal = "IID"
)
m2

tidy(m2, conf.int = TRUE)
tidy(m2, "ran_pars", conf.int = TRUE)

## tweedie density model ------

m3 <- sdmTMB(
  data = pcod,
  formula = density ~ poly(log(depth), 2),
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on",
  time = "year",
  spatiotemporal = "IID"
)
m3

tidy(m3, conf.int = TRUE)
tidy(m3, "ran_pars", conf.int = TRUE)  # Tweedie p (power) parameter; between 1 and 2.


   #  AR1- influence of past decays over time and process is stationary in that it tends to revert to the mean

   # takes a little bit longer.... about a minute on my laptop
fit_spatiotemporal <- sdmTMB(
  density ~ s(depth, k = 5), 
  data = pcod, 
  mesh = mesh,
  time = "year",
  family = tweedie(link = "log"), 
  spatial = "off", 
  spatiotemporal = "ar1"   # can also be 'iid' or 'rw' (or 'off')
)

table(pcod$year)
fit_spatiotemporal

tidy(fit_spatiotemporal, conf.int = TRUE)
tidy(fit_spatiotemporal, "ran_pars", conf.int = TRUE)  # Tweedie p (power) parameter; between 1 and 2.

pcod$resids <- residuals(fit_spatiotemporal, type = "mle-mvn") # randomized quantile residuals
qqnorm(pcod$resids)
abline(0, 1)

ggplot(pcod, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~year) +
  coord_fixed()

s <- simulate(fit_spatiotemporal, nsim = 300, type = "mle-mvn")
dharma_residuals(s, fit_spatiotemporal)


## index of predicted abundance/density over time -----

grid_yrs <- replicate_df(qcs_grid, "year", unique(pcod$year))
p_st <- predict(fit_spatiotemporal, newdata = grid_yrs, 
                return_tmb_object = TRUE)
index <- get_index(p_st, area = rep(4, nrow(grid_yrs)))
ggplot(index, aes(year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
  geom_line(lwd = 1, colour = "grey30") +
  labs(x = "Year", y = "Biomass (kg)") +
  theme_classic()


## center of gravity (for detecting distribution shifts...)  --------

cog <- get_cog(p_st, format = "wide")
ggplot(cog, aes(est_x, est_y, colour = year)) +
  geom_pointrange(aes(xmin = lwr_x, xmax = upr_x)) +
  geom_pointrange(aes(ymin = lwr_y, ymax = upr_y)) +
  scale_colour_viridis_c()


# time varying cofficient --------

# NOTE: with 'rw' (default) you should not specify the same params in the main and time varying formulas

fit <- sdmTMB(
  density ~ 0 + s(depth, k = 5), 
  time_varying = ~ 1,      # time varying intercept
  data = pcod, 
  mesh = mesh,
  time = "year", # NOTE: this automatically makes this a spatiotemporal model 
  family = tweedie(link = "log"),
  silent = FALSE, # see progress
  spatial="off"
)

fit

# ?tidy.sdmTMB
tidy(fit,effects='fixed')
tidy(fit,effects='ran_pars')
tidy(fit,effects='ran_vals')   # random intercept

grid_yrs <- replicate_df(qcs_grid, time_name = "year", time_values = unique(pcod$year))
p <- predict(fit, newdata = grid_yrs)

names(p)
nrow(p) / nrow(qcs_grid)   # note that this dataset replicates qcs_grid 9 times (9 years)

ggplot(p,aes(X,Y,fill=est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)

names(p)
ggplot(p,aes(X,Y,fill=est_rf)) +  # all random fields compbined
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)
  
names(p)  # note no omega: only fits stf if iid
ggplot(p,aes(X,Y,fill=epsilon_st)) +  # spatiotemporal field: Note different spatial field for each year
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)

names(p)
ggplot(p,aes(X,Y,fill=est_non_rf)) +  # spatial variation in predictions driven by depth (not temporally varying)
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)

## where does the time varying intercept come in?

tapply(p$est_non_rf,p$year,mean)  # no difference in non-rf predictions by year  (year effect is a random effect)

tapply(p$epsilon_st,p$year,mean)   # I guess this is where you can 'see' the random time-varying intercept...

## time-varying effect of depth -------

   ## note: can't include a smooth in time-varying coefs
          # can't have the same variables in main formula and time-varying formula- including intercept

fit <- sdmTMB(
  density ~ 1,   # global intercept
  time_varying = ~ 0 + depth_scaled + depth_scaled2,  # time varying effect of depth
  data = pcod, mesh = mesh,
  time = "year",
  family = tweedie(link = "log"),
  spatial = "off",
  spatiotemporal = "ar1",
  silent = FALSE
)

grid_yrs <- replicate_df(qcs_grid, time_name = "year", time_values = unique(pcod$year))
p <- predict(fit, newdata = grid_yrs)

names(p)
ggplot(p,aes(X,Y,fill=est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)

names(p)
ggplot(p,aes(X,Y,fill=est_rf)) +  # all random fields compbined
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)

names(p)  
ggplot(p,aes(X,Y,fill=epsilon_st)) +  # spatiotemporal field: Note different spatial field for each year
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)

names(p)
ggplot(p,aes(X,Y,fill=est_non_rf)) +  # spatial variation in predictions driven by depth (now time-varying)
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)


## try AR1 time varying model with spatial field only

fit <- sdmTMB(
  density ~ 1 + depth_scaled + depth_scaled2,   # global intercept
  time_varying = ~ 1 + depth_scaled + depth_scaled2,  # time varying effect of depth
  time_varying_type = "ar1",
  data = pcod, mesh = mesh,
  extra_time = seq(min(pcod$year), max(pcod$year)),
  family = tweedie(link = "log"),
  spatial = "on",
  time = "year",
  spatiotemporal = "off",
  silent=FALSE
)

fit

nd <- expand.grid(
  depth_scaled = seq(min(pcod$depth_scaled) + 0.2,
                     max(pcod$depth_scaled) - 0.2,
                     length.out = 50
  ),
  year = unique(pcod$year) # all years
)
nd$depth_scaled2 <- nd$depth_scaled^2

p <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA)

ggplot(p, aes(depth_scaled, exp(est),
              ymin = exp(est - 1.96 * est_se),
              ymax = exp(est + 1.96 * est_se),
              group = as.factor(year)
)) +
  geom_line(aes(colour = year), lwd = 1) +
  geom_ribbon(aes(fill = year), alpha = 0.1) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  scale_x_continuous(labels = function(x) round(exp(x * pcod$depth_sd[1] + pcod$depth_mean[1]))) +
  coord_cartesian(expand = F) +
  labs(x = "Depth (m)", y = "Biomass density (kg/km2)")

# spatially varying coefficient ------------

## time trend varies by space ------

pcod$year_scaled <- as.numeric(scale(pcod$year))
fit <- sdmTMB(
  density ~ s(depth, k = 5) + year_scaled,
  spatial_varying = ~ year_scaled,   # note that includes both global effect of year and spatial-varying effect
  data = pcod, 
  mesh = mesh, 
  time = "year",   # is this needed?
  family = tweedie(link = "log"),
  spatiotemporal = "off",  
  spatial="on"
)

fit 

allyears=sort(unique(pcod$year)); allyears

grid_yrs <- replicate_df(qcs_grid, "year", unique(pcod$year))
grid_yrs$year_scaled <- (grid_yrs$year - mean(pcod$year)) / sd(pcod$year)
p <- predict(fit, newdata = grid_yrs) %>% 
  subset(year == 2003) # any year  - doesn't matter...
ggplot(p, aes(X, Y, fill = zeta_s_year_scaled)) + geom_raster() +
  scale_fill_gradient2()


# random intercepts ----

  ## note that random slopes a la lme4 are not allowed. But can fit factor smoothers and 2-d smoothers

pcod$year_factor <- as.factor(pcod$year)
fit <- sdmTMB(
  density ~ s(depth, k = 5) + (1 | year_factor),   # random intercept for year
  data = pcod, mesh = mesh,
  time = "year",
  family = tweedie(link = "log")
)

fit  # note: spatiotemporal model, random intercept for year (same as time varying intercept?)


# breakpoints and thresholds -----------

fit <- sdmTMB(
  present ~ 1 + breakpt(depth_scaled), 
  data = pcod, mesh = mesh,
  family = binomial(link = "logit")
)

fit <- sdmTMB(
  present ~ 1 + logistic(depth_scaled), 
  data = pcod, mesh = mesh,
  family = binomial(link = "logit")
)


# simulate data from scratch ---------

predictor_dat <- expand.grid(
  X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100)
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.05)
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = poisson(link = "log"),
  range = 0.3,
  sigma_O = 0.4,
  seed = 1,
  B = 1 # B0 = intercept
)
head(sim_dat)

  # sample 200 points for fitting
sim_dat_obs <- sim_dat[sample(seq_len(nrow(sim_dat)), 200), ]


ggplot(sim_dat, aes(X, Y)) +
  geom_raster(aes(fill = exp(eta))) + # mean without observation error
  geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)

## fit to simulated data -------

mesh <- make_mesh(sim_dat_obs, xy_cols = c("X", "Y"), cutoff = 0.05)
fit <- sdmTMB(
  observed ~ 1,
  data = sim_dat_obs,
  mesh = mesh,
  family = poisson()
)

fit
p <- predict(fit,newdata = sim_dat)
head(p)
names(p)

ggplot(p,aes(X,Y,fill=est))+geom_raster()+scale_fill_viridis_c()


# simulate from existing fit ----------

s <- simulate(fit, nsim = 500)
dim(s)  # 500 sims of 200 observed data

s[1:3,1:4]


# uncertainty in parameters and derived params ----------

samps <- gather_sims(fit, nsim = 1000)
dim(samps)   # 1 sample for each of three parameters
names(samps)
head(samps)

unique(samps$.variable)

ggplot(samps, aes(.value)) + geom_histogram() +
  facet_wrap(~.variable, scales = "free_x")


# uncertainty on spatial predictions ----------

p <- predict(fit, newdata = predictor_dat, nsim = 500)
dim(p)
predictor_dat$se <- apply(p, 1, sd)  # uncertainty 
ggplot(predictor_dat, aes(X, Y)) +
  geom_raster(aes(fill = se)) +
  scale_fill_viridis_c(option = "A") +
  coord_cartesian(expand = FALSE) +
  geom_point(data=sim_dat_obs,aes(X,Y))

plot(mesh)


# cross validation ----------

data(pcod)
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
pcod$fyear <- as.factor(pcod$year)

plot(mesh)

## Set parallel processing if desired:
library(future)
plan(multisession, workers = 2)

m_cv <- sdmTMB::sdmTMB_cv(
  density ~ 0 + s(depth_scaled) + fyear,
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  k_folds = 4
) 
closeAllConnections()

m <- sdmTMB(
  density ~ 0 + s(depth_scaled) + fyear,
  data = pcod, mesh = mesh,
  family = tweedie(link = "log")
)
m


# names(m_cv$data)  # add columns for cv_predicted, etc.
p = predict(m)

# Sum of log likelihoods of left-out data: equivalent to ELPD - useful for model comparison!
m_cv$sum_loglik   # compare with -6403.637 for fitted model

m_cv$fold_loglik  # log likelihood of each fold

# RMSE across entire dataset:
sqrt(mean((pcod$density - exp(p$est))^2))    # fitted
sqrt(mean((pcod$density - m_cv$data$cv_predicted)^2))   # CV
sqrt(mean((m_cv$data$density - m_cv$data$cv_predicted)^2))    # CV alternative

# MAE across entire dataset:
mean(abs(pcod$density-exp(p$est)))   # fitted
mean(abs(m_cv$data$density - m_cv$data$cv_predicted))   # CV

### ROC 

library(pROC)
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(present ~ s(depth), data = pcod, mesh = mesh)
pred <- predict(fit) # presence-absence model
roc <- pROC::roc(pcod$present, plogis(pred$est))
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases
auc <- pROC::auc(roc)
auc


x <- sdmTMB_cv(
  present ~ s(depth), data = pcod, spatial = "off",
  mesh = mesh, family = binomial(), k_folds = 2
)

roc <- pROC::roc(x$data$present, plogis(x$data$cv_predicted))
auc <- pROC::auc(roc)
auc


## spatial clustering -------
clust <- kmeans(pcod[, c("X", "Y")], 20)$cluster

m_cv <- sdmTMB::sdmTMB_cv(
  density ~ 0 + s(depth_scaled) + fyear,
  data = pcod,
  mesh = mesh,
  fold_ids = clust,
  family = tweedie(link = "log")
)

# priors ------------

mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  density ~ depth_scaled,
  data = pcod, mesh = mesh,
  family = tweedie(),
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 10, sigma_lt = 5),
    b = normal(c(0, 0), c(1, 10)),
    phi = halfnormal(0, 15)
  )
)

plot_pc_matern(range_gt = 10, sigma_lt = 5)


# custom mesh ----------

bnd <- fmesher::fm_nonconvex_hull(cbind(pcod$X, pcod$Y), convex = -0.1)
mesh_inla <- fmesher::fm_mesh_2d_inla(
  boundary = bnd,
  max.edge = c(25, 50)
)
mesh <- make_mesh(pcod, c("X", "Y"), mesh = mesh_inla)
plot(mesh)


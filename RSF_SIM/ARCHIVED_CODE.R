
## SIMPLIFIED - ZINB with Spatial Smooth in gamlss
## Essential code for quick implementation

library(gamlss)
library(gamlss.add)
library(DHARMa)

## Basic ZINB model with 2D spatial smooth -------

# For count data (negative binomial)
crest6 <- gamlss(
  cbind(lark, nsurv - lark) ~ ga(~s(east, north, k = 100)),  # 2D P-spline spatial smooth
  sigma.formula = ~ 1,               # constant dispersion
  nu.formula = ~ 1,                  # constant zero-inflation
  family = ZIBI,                    # Zero-Inflated Negative Binomial
  data = bird2
)

summary(crest6)

## Quick diagnostics -------
plot(crest6)  # standard gamlss plots
wp(zibi_model)     # worm plot

## DHARMa -------
sim_matrix <- matrix(NA, nrow = nrow(bird2), ncol = 1000)
for(i in 1:1000) {
  sim_matrix[, i] <- rZINBI(
    n = nrow(bird2),
    mu = fitted(crest6, what = "mu"),
    sigma = fitted(crest6, what = "sigma")
    # nu = fitted(crest6, what = "nu")
  )
}

res <- createDHARMa(
  simulatedResponse = sim_matrix,
  observedResponse = bird2$lark,
  fittedPredictedResponse = fitted(crest6),
  integerResponse = TRUE
)

plot(res)
testZeroInflation(res)
testSpatialAutocorrelation(res, x = bird2$east, y = bird2$north)

## Predictions -------
new_data <- data.frame(east = 50, north = 50)
predict(zinb_model, newdata = new_data, type = "response")

## Key Parameters -------

# pb(east, north, df = X) - 2D spatial smooth
#   df: degrees of freedom (higher = more flexible)
#   typical values: 9, 16, 25

# Families:
#   ZINBI - Zero-Inflated Negative Binomial (for counts)
#   ZIBI  - Zero-Inflated Binomial (for proportions)
#   ZIP   - Zero-Inflated Poisson
#   NBI   - Negative Binomial (no ZI)

# Formulas:
#   mu.formula (main) - mean structure
#   sigma.formula     - dispersion parameter
#   nu.formula        - zero-inflation probability
#   tau.formula       - fourth parameter (if available)

## Alternative smooths -------
# pb() - P-splines (recommended for spatial)
# lo() - loess smooth
# cs() - cubic spline
# nn() - neural network




library(sdmTMB)

# Scale coordinates (recommended for GP models)
bird2$east_scaled <- scale(bird2$east)[,1]
bird2$north_scaled <- scale(bird2$north)[,1]

mesh <- make_mesh(bird2, xy_cols = c("east_scaled", "north_scaled"), cutoff = 0.1)

plot(mesh)

fit <- sdmTMB(
  cbind(lark, nsurv - lark) ~ 1,
  data = bird2,
  mesh = mesh,
  family = binomial()
  # spatial = "on"
)

summary(fit)
sanity(fit)
plot(fit)








## SIMPLIFIED VERSION - Zero-Inflated Binomial with GP in brms
## Essential code for quick reference

library(brms)
library(DHARMa)

## Prepare data -------



## Fit the model -------

## simple brms model

crest1 <- brm(lark | trials(nsurv) ~ 1,
              cores=4,backend="cmdstanr",
              data = bird2, family = binomial() )

summary(crest2)

crest_brms1 <- brm(lark | trials(nsurv) ~ 1,
                   cores=4,backend="cmdstanr",
                   data = bird2, family = binomial() )

crest_brms2 <- brm(
  lark | trials(nsurv) ~ 1 + gp(east_scaled, north_scaled,k=100),
  cores=4,backend="cmdstanr",
  data = bird2, family = binomial() 
)

crest_brms3 <- brm(
  lark | trials(nsurv) ~ 1 + gp(east_scaled, north_scaled,k=100),
  cores=4,backend="cmdstanr",
  data = bird2, family = zero_inflated_binomial()  
)

#   # Zero-inflation component
#   zi ~ 1,
#   
#   # Family
#   family = zero_inflated_binomial(),
#   
#   # Data
#   data = bird2,
#   
#   # MCMC settings
#   chains = 4,
#   iter = 2000,
#   warmup = 1000,
#   cores = 4,
#   backend = "cmdstanr",
#   
#   # Help with convergence
#   control = list(adapt_delta = 0.95),
#   
#   # Cache the model
#   file = "crest_brms_model"
# )

## Check the results -------
summary(crest_brms)

# Check convergence
rhat(crest_brms)  # should be < 1.01
neff_ratio(crest_brms)  # should be > 0.1

# Posterior predictive check
pp_check(crest_brms, ndraws = 100)

## DHARMa diagnostics -------
res <- createDHARMa(
  simulatedResponse = t(posterior_predict(crest_brms, ndraws = 1000)),
  observedResponse = bird2$lark,
  fittedPredictedResponse = apply(
    t(posterior_predict(crest_brms, ndraws = 1000)), 1, median
  ),
  integerResponse = TRUE
)

plot(res)
testZeroInflation(res)
testSpatialAutocorrelation(res, x = bird2$east, y = bird2$north)

## Make predictions -------
# Predict for new data
new_data <- data.frame(
  east_scaled = 0,
  north_scaled = 0,
  nsurv = 1
)

# Get posterior predictions
posterior_epred(crest_brms, newdata = new_data)

## Notes -------
# - gp(east_scaled, north_scaled) creates Gaussian process spatial effect
# - zi ~ 1 adds zero-inflation with just an intercept
# - Use scaled coordinates for better MCMC sampling
# - file = "name" caches the model to avoid refitting
# - Increase adapt_delta if you get divergent transitions
# - Use posterior_epred() for expected predictions
# - Use posterior_predict() for predicted observations (with noise)














library(glmmTMB)
library(sdmTMB)
library(ropenblas)

m <- 1e4; n <- 1e3; k <- 3e2
X <- matrix(rnorm(m*k), nrow=m); Y <- matrix(rnorm(n*k), ncol=n)
system.time(X %*% Y)

library(dplyr)
library(ggplot2)
library(sdmTMB)
data(pcod)
head(pcod)

mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)

plot(mesh)

fit <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on"
)

summary(fit)

tidy(fit, conf.int = TRUE)
tidy(fit, effects = "ran_pars", conf.int = TRUE)
sanity(fit)

ggeffects::ggpredict(fit, "depth [50:400, by=2]") |> plot()

p <- predict(fit, newdata = qcs_grid)

ggplot(p, aes(X, Y, fill = exp(est))) + geom_raster() +
  scale_fill_viridis_c(trans = "sqrt")




# n <- 25                                              ## Number of time points
# x <- MASS::mvrnorm(mu = rep(0,n),
#                    Sigma = .7 ^ as.matrix(dist(1:n)) )    ## Simulate the process using the MASS package
# y <- x + rnorm(n)
# times <- factor(1:n, levels=1:n)
# head(levels(times))
# group <- factor(rep(1,n))
# dat0 <- data.frame(y, times, group)
# glmmTMB(y ~ ar1(times + 0 | group), data=dat0)

# Scale coordinates for better model convergence (optional)

bird2$east_scaled <- scale(bird2$east)[,1]
bird2$north_scaled <- scale(bird2$north)[,1]
bird2$pos_scaled <- numFactor(bird2$east_scaled, bird2$north_scaled)
bird2$group = 1
bird2$lark
table(bird2$nsurv)

## Fit zero-inflated binomial model with Gaussian process spatial effect -------
# Using scaled coordinates

crest_zi <- glmmTMB(
  cbind(lark, nsurv - lark) ~ 1 + mat(pos_scaled+0|group),
  ziformula = ~ 1,  # zero-inflation intercept only
  family = binomial(),
  data = bird2
)

summary(crest_zi)


## Model checking with DHARMa -----
cat("\n=== DHARMa Diagnostics for Zero-Inflated Model ===\n")

# Simulate residuals
res_zi <- simulateResiduals(crest_zi, n = 1000)

# Overall diagnostic plot
plot(res_zi, main = "Zero-Inflated Binomial with GP Spatial Effect")

# Test residuals
testResiduals(res_zi)

# Test for remaining zero-inflation
cat("\n--- Test for Zero-Inflation ---\n")
test_zi <- testZeroInflation(res_zi)
print(test_zi)

# Test for outliers
cat("\n--- Test for Outliers ---\n")
test_outliers <- testOutliers(res_zi, type = "bootstrap")
print(test_outliers)

# Test for dispersion
cat("\n--- Test for Dispersion ---\n")
test_disp <- testDispersion(res_zi)
print(test_disp)

# Test for spatial autocorrelation
cat("\n--- Test for Spatial Autocorrelation ---\n")
test_spatial <- testSpatialAutocorrelation(
  res_zi,
  x = bird2$east,
  y = bird2$north
)
print(test_spatial)

# Test for uniformity
cat("\n--- Test for Uniformity ---\n")
test_unif <- testUniformity(res_zi)
print(test_unif)

## Additional diagnostic plots -----
# Residuals vs fitted
ggplot(data.frame(Fitted = fitted(crest_zi),
                  Resid = residuals(res_zi)),
       aes(Fitted, Resid)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "DHARMa Residuals vs Fitted Values",
       x = "Fitted Values", y = "DHARMa Residuals")

# Spatial pattern of residuals
bird2$dharma_resid <- residuals(res_zi)

ggplot(bird2, aes(x = east, y = north, color = dharma_resid)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0.5,
    name = "DHARMa\nResidual"
  ) +
  theme_minimal() +
  labs(title = "Spatial Pattern of DHARMa Residuals",
       subtitle = "Zero-Inflated Binomial with GP Spatial Effect",
       x = "East", y = "North")

## Model Comparison -------
cat("\n=== Model Comparison ===\n")

# 1. Standard binomial with GP (no zero-inflation)
crest_binom <- glmmTMB(
  cbind(lark, nsurv - lark) ~ 1 + gp(pos_scaled),
  family = binomial,
  data = bird2
)

# 2. Zero-inflated binomial without spatial effect
crest_zi_nospatial <- glmmTMB(
  cbind(lark, nsurv - lark) ~ 1,
  ziformula = ~ 1,
  family = binomial,
  data = bird2
)

# 3. Standard binomial without spatial effect (simplest model)
crest_simple <- glmmTMB(
  cbind(lark, nsurv - lark) ~ 1,
  family = binomial,
  data = bird2
)

# AIC comparison
aic_table <- AIC(crest_zi, crest_binom, crest_zi_nospatial, crest_simple)
aic_table <- aic_table[order(aic_table$AIC),]
cat("\nAIC Comparison (lower is better):\n")
print(aic_table)

# DHARMa diagnostics for binomial model without ZI
cat("\n=== DHARMa Diagnostics for Standard Binomial (no ZI) ===\n")
res_binom <- simulateResiduals(crest_binom, n = 1000)
plot(res_binom, main = "Binomial (no ZI) with GP Spatial Effect")
testZeroInflation(res_binom)

# DHARMa diagnostics for ZI model without spatial effect
cat("\n=== DHARMa Diagnostics for ZI Model (no spatial effect) ===\n")
res_zi_nospatial <- simulateResiduals(crest_zi_nospatial, n = 1000)
plot(res_zi_nospatial, main = "Zero-Inflated Binomial (no spatial effect)")
testSpatialAutocorrelation(res_zi_nospatial, x = bird2$east, y = bird2$north)

## Extract model parameters -------
cat("\n=== Model Parameters (Best Model) ===\n")

# Fixed effects
cat("\nFixed Effects:\n")
print(fixef(crest_zi))

# Confidence intervals
cat("\nConfidence Intervals:\n")
print(confint(crest_zi))

# Variance of Gaussian process
cat("\nVariance Components:\n")
print(VarCorr(crest_zi))

# Zero-inflation probability
zi_prob <- plogis(fixef(crest_zi)$zi)
cat("\nEstimated zero-inflation probability:", round(zi_prob, 3), "\n")

## Predictions and visualization -------
# Create prediction grid
east_seq <- seq(min(bird2$east), max(bird2$east), length.out = 50)
north_seq <- seq(min(bird2$north), max(bird2$north), length.out = 50)
pred_grid <- expand.grid(east = east_seq, north = north_seq)

# Need to create scaled versions and pos for predictions
pred_grid$east_scaled <- (pred_grid$east - mean(bird2$east)) / sd(bird2$east)
pred_grid$north_scaled <- (pred_grid$north - mean(bird2$north)) / sd(bird2$north)
pred_grid$pos_scaled <- numFactor(pred_grid$east_scaled, pred_grid$north_scaled)
pred_grid$nsurv <- 1  # predictions for single survey

# Get predictions (on probability scale)
pred_grid$pred <- predict(crest_zi, newdata = pred_grid, type = "response")

# Plot predicted probability surface
ggplot(pred_grid, aes(x = east, y = north, fill = pred)) +
  geom_tile() +
  geom_point(data = bird2, aes(x = east, y = north), 
             inherit.aes = FALSE, alpha = 0.3, size = 1) +
  scale_fill_viridis_c(name = "Predicted\nProbability") +
  theme_minimal() +
  labs(title = "Predicted Probability of Crested Lark Presence",
       subtitle = "Black points show observed locations",
       x = "East", y = "North")

cat("\n=== Analysis Complete ===\n")


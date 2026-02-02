# Reanalyse arthropod data from Seibold et al (2019) using some
# ideas from Daskalova, Phillimore & Meyers-Smith (2021)
#
# Refs
# Daskalova, G.N., Phillimore, A.B., Myers‐Smith, I.H., 2021. Accounting for
#   year effects and sampling error in temporal analyses of invertebrate
#   population and biodiversity change: a comment on Seibold et al. 2019.
#   *Insect Conserv. Divers.* 14, 149–154. https://doi.org/10.1111/icad.12468

# Seibold, S., Gossner, M.M., Simons, N.K., Blüthgen, N., Müller, J., Ambarlı,
#   D., Ammer, C., Bauhus, J., Fischer, M., Habel, J.C., Linsenmair, K.E.,
#   Nauss, T., Penone, C., Prati, D., Schall, P., Schulze, E.-D., Vogt, J.,
#   Wöllauer, S., Weisser, W.W., 2019. Arthropod decline in grasslands and
#   forests is associated with landscape-level drivers. Nature 574, 671–674.
#   https://doi.org/10.1038/s41586-019-1684-3

rm(list=ls())

# load packages ---------

pkgs <- c("here", "readr", "janitor", "mgcv", "gratia", "dplyr", "ggplot2",
          "ggrepel")
vapply(pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)

# read data -------

  # data - from Seibold et al Nature 2019
  # https://doi.org/10.25829/bexis.25786-1.3.11  NOTE: column names are slightly different in this download
seibold <- read_csv2(here("data/25786_3_data.csv"),
                     col_types = "cccccdndddddcddnddcddcddcddnnnnnnnnnnnn")

## basic data wrangling -------

seibold <- seibold %>%
  janitor::clean_names() %>%   
  rename(year = collection_year,
         region = exploratory,
         habitat = habitat_type) %>%
  mutate(across(c(habitat, sampling_regime, region, plot_id_year, plot_id),
                as.factor))

# visualization and exploration -------

seibold %>%
  ggplot(aes(x = year, y = abundance_identified, group = plot_id)) +
  geom_line() +
  facet_grid(region ~ habitat)

## more wrangling
seibold <- seibold %>%
  mutate(year_f = factor(year), # year as a factor for ranef
         year_c = year - 2012) # centre year

# how many plots
seibold %>%
  group_by(plot_id) %>%
  count()

seibold %>%
  group_by(habitat, region, plot_id) %>%
  count() %>%
  arrange(-n)

# filter out the forest sites; many have fewer than 9 observations
seibold <- seibold %>%
  filter(habitat == "grassland")

seibold_end <- seibold %>%
  group_by(plot_id) %>%
  top_n(1, year) %>%
  ungroup() %>%
  mutate(year = max(year))

## data subset for GAM test --------

seibold %>%
  filter(plot_id == "AEG1") %>%
  ggplot(aes(x = year, y = abundance_identified)) +
  geom_line(alpha = 0.5) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 9)) +
  labs(y = "Abundance", x = NULL)

site <- seibold %>%
  filter(plot_id == "AEG1")


# GAM model ------------


   # Negative binomial with theta (constant) estimated
m_site <- gam(abundance_identified ~ s(year), data = site,
              method = "REML", family = nb())

summary(m_site) # model summary

draw(m_site) # partial effect of smooth (use gratia package for nice plot)

## model checking --------

# plot observed vs predicted
plot(site$abundance_identified, predict(m_site,type="response"))
abline(1,1)

# fancy method in gratia package! 

appraise(m_site, method = "simulate") # visual model diagnostics (gratia)

k.check(m_site) # check basis size (p value should be relatively high for all basis expansions)

# plot residuals over time (always good to look for patterns in the residuals!)
site %>%
  add_residuals(m_site) %>%    # gratia function
  ggplot(aes(x = year, y = .residual)) +
    geom_line() +
    geom_point()

acf(residuals(m_site))  # check for temporal patterns in residuals

library(DHARMa)
res = simulateResiduals(m_site)
testResiduals(res)

# fit model to full dataset -------

# use multiple threads for fitting
ctrl <- gam.control(nthreads = 4)

seibold$region

# fit to all data (takes a couple minutes)
m1 <- gam(abundance_identified ~ s(year,k=5) + # smooth of year (trend)
            s(region,year,bs="sz",k=5,id=1) +  # allow year effect to vary by region
            s(plot_id, bs = "re"), # random effect of plot
          data = seibold,
          method = "REML",
          family = nb(),
          control = ctrl)

?smooth.construct.sz.smooth.spec

res1=DHARMa::simulateResiduals(m1)
testResiduals(res1)


plot(m1)

draw(m1) 

summary(m1)

library(ggeffects)
pred = ggpredict(m1, terms = "year [all]")
plot(pred)

summary(m1)
pred <- ggpredict(m1, terms = c("year [all]", "region [all]"))
plot(pred)

names(sample(which(table(seibold$plot_id)>0),5))
pred <- ggpredict(m1,terms = c("year [all]", "plot_id [SEG32,HEG42,SEG10]") )
plot(pred)
  

# alternative formula for region-specific trends  - also using tweedie error distribution
m2 <- gam(abundance_identified ~ region + # regional means
            s(year, by = region) + # region-specific trends (no main effect)
            s(plot_id, bs = "re"), # random intercept for plot
          data = seibold,
          method = "REML",
          family = tw(),
          control = ctrl)

summary(m2)
draw(m2)

# individual trends?
# about 1 minute 
# system.time(
# m3 <- gam(abundance_identified ~
#             s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
#           data = seibold,
#           method = "REML",
#           family = nb(),
#           control = ctrl)
#   )

# summary(m3)
# 
# draw(m3)


# region-specific trends, plus individual trends, plus year-to-year effects
# about 5 minutes on my laptop
# system.time(
# m5 <- gam(abundance_identified ~ region + # regional means fixef
#             s(year_f, bs = "re") + # year-to-year effects
#             s(year, by = region, k = 20) + # region-specific smooths
#             s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
#           data = seibold,
#           method = "REML",
#           family = nb(),
#           control = ctrl)
#   )

# # plot smooths
# draw(m5)
# 
# # model diagnostics
# appraise(m5, method = "simulate")

# basis size
# k.check(m5)

# rootogram
# rootogram(m5, max_count = 300) %>%
  # draw()

## Alternatives with bam()
## not threaded unless you set up a cluster -- see ?bam
## DON'T RUN THIS IN WORKSHOP!! - about 17 minutes
# system.time(
# b5a <- bam(abundance_identified ~ region + # regional means fixef
#              s(year_f, bs = "re") + # year-to-year effects
#              s(year, by = region) + # region-specific smooths
#              s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
#            data = seibold,
#            method = "fREML", # <-- fast REML!
#            family = nb(),
#            control = ctrl)
#   )

## Discretise covariates -- algorithm uses `nthreads` threads
# system.time(
# b5b <- bam(abundance_identified ~ region + # regional means fixef
#              s(year_f, bs = "re") + # year-to-year effects
#              s(year, by = region) + # region-specific smooths
#              s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
#            data = seibold,
#            method = "fREML", # <-- fast REML!
#            discrete = TRUE,  # <-- discretise covariates == smaller basis
#            family = nb(),
#            control = ctrl)
#   )

## Location scale models
## about 10 minutes!!
# system.time(
# m_twlss <- gam(list(
#                  abundance_identified ~ region + # regional means fixef
#                    s(year_f, bs = "re") + # year-to-year effects
#                    s(year, by = region) + # region-specific smooths
#                    s(year, plot_id, bs = "fs", k = 5), # plot-specific trends
#                  ~ s(plot_id, bs = "re"),
#                  ~ s(plot_id, bs = "re")
#                 ),
#               data = seibold,
#               method = "REML",
#               optimizer = "efs",
#               family = twlss(),
#               control = ctrl)
#   )

# summary(m_twlss)
# 
# draw(m_twlss)
# 
# appraise(m_twlss)

##
# smooths(m5)
# y2y <- smooth_estimates(m5, "s(year_f)")
# 
# y2y %>%
#   ggplot(aes(x = year_f, y = est)) +
#   geom_pointrange(aes(ymin = est - se,
#                       ymax = est + se)) +
#   labs(x = NULL)

## other covariates
# seibold <- seibold %>%
#   mutate(year_s = scale(year),
#          landuse_intensity_s = scale(landuse_intensity)[,1],
#          mean_winter_temperature_s = scale(mean_winter_temperature)[,1],
#          precipitation_sum_growing_preriod_s =
#            scale(precipitation_sum_growing_preriod)[,1],
#          grassland_cover_1000_s = scale(grassland_cover_1000)[,1],
#          arable_cover_1000_s = scale(arable_cover_1000)[,1])

## about 1.5 minutes
# system.time(
# m_cov <- gam(abundance_identified ~
#                s(year_s) + # overall trend
#                s(landuse_intensity_s) +
#                s(mean_winter_temperature_s) +
#                s(precipitation_sum_growing_preriod_s) +
#                s(grassland_cover_1000_s) +
#                s(arable_cover_1000_s) +
#                ti(mean_winter_temperature_s,
#                   precipitation_sum_growing_preriod_s) +
#                ti(year_s, landuse_intensity_s) +
#                ti(year_s, grassland_cover_1000_s) +
#                ti(year_s, arable_cover_1000_s) +
#                s(year_f, bs = "re") + # year-to-year effects
#                s(plot_id, bs = "re"), # site specific mean abundance
#              family = nb(),
#              method = "REML",
#              control = ctrl,
#              data = seibold,
#              select = TRUE)
#   )

# plot the smooths
# sms <- smooths(m_cov)

# plot the univariate smooths
# draw(m_cov, select = sms[1:6])

# plot the tensor product interation smooths
# draw(m_cov, select = sms[7:10], rug = FALSE)

# plot the ranefs
# draw(m_cov, select = sms[11:12])

# plot the overall trend effect
# draw(m_cov, select = "s(year_s)")

# DON'T RUN IN WORKSHOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# put in site specific trends
# ~17 minutes 
# system.time(
# m_cov2 <- gam(abundance_identified ~
#                 s(year_s) +
#                 s(landuse_intensity_s) +
#                 s(mean_winter_temperature_s) +
#                 s(precipitation_sum_growing_preriod_s) +
#                 s(grassland_cover_1000_s) +
#                 s(arable_cover_1000_s) +
#                 ti(mean_winter_temperature_s,
#                    precipitation_sum_growing_preriod_s) +
#                 ti(year_s, landuse_intensity_s) +
#                 ti(year_s, grassland_cover_1000_s) +
#                 ti(year_s, arable_cover_1000_s) +
#                 s(year_f, bs = "re") +
#                 s(year_s, plot_id, bs = "fs", k = 5), # <-- here
#               family = nb(),
#               method = "REML",
#               control = ctrl,
#               data = seibold)
#   )

# AIC(m_cov, m_cov2)

##
# smooths(m_cov2)
# y2y <- smooth_estimates(m_cov2, "s(year_f)")

# y2y %>%
#   ggplot(aes(x = year_f, y = est)) +
#   geom_pointrange(aes(ymin = est - se,
#                       ymax = est + se)) +
#   labs(x = NULL)
# 
# draw(m_cov2, select = "s(year_s,plot_id)")

# year, day-of-year

# knots <- list(day_of_year = c(0.5, 366.5))
# gam(y ~ s(year) + s(day_of_year, bs = "cc"), knots = knots)

# knots <- list(month = c(0.5, 12.5))
# gam(y ~ s(year) + s(month, bs = "cc", k = 12), knots = knots)

# gam(y ~ te(year, month, bs = c("tp", "cc")), knots = knots)
# gam(y ~ s(year) + s(month, bs = "cc", k = 12) +
#       ti(year, month, bs = c("tp", "cc")), knots = knots)


# END SCRIPT


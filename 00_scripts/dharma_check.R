# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

library(DHARMa)

# due to our large sample size, tests are bound to turn up significant
# we need to look at magnitude of change and decide whether or not it is a problem



# setup
anals <- c("d02_fidelity_po", "d03_time_po", "d06_distance", "d07_duration", 
           "d08_length",
           "d01_group_po",  "d04_hotspot", "d05_protocol", "d09b_s_cover", 
           "d09a_s_UNU")

anal_name <- anals[10]

load(glue("outputs/{anal_name}.RData"))




# simulated residuals
tictoc::tic("Simulating residuals for model")
sim_res <- simulateResiduals(fittedModel = model_a, plot = F)
# sim_res2 <- simulateResiduals(fittedModel = model_b, plot = F)
tictoc::toc()


# plotting residuals
plot(sim_res)
plot(sim_res2)


# testing
testUniformity(sim_res)
# testOutliers(sim_res)
testDispersion(sim_res)
# testTemporalAutocorrelation(recalculateResiduals(sim_res, group = data_a$MONTH), 
#                             time = unique(data_a$MONTH))

# testing for heteroscedasticity
plotResiduals(sim_res, data_a$MONTH)
plotResiduals(sim_res, data_a$M.YEAR)

# plotResiduals(sim_res, (data_a %>% filter(!is.na(COUNTY)))$MONTH)
# plotResiduals(sim_res, (data_a %>% filter(!is.na(COUNTY)))$M.YEAR)




# Most models: some variation of the middle region of the QQ is either higher or lower 
# than the expected values. According to DHARMa vignette this is not the usual 
# over/underdispersion pattern, but rather suggest heteroscedasticity. 
# Probably because we only have month and year in our models and no other variable 
# that could probably have a big effect on the response variables.
# plot of response vs predictors: individual levels of the predictor are pretty similar, 
# and just their means are either slightly higher or lower than the dotted line.
# Should be okay to ignore.
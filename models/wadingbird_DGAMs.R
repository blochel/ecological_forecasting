
# DGAM models - Everglades Wading Birds -----------------------------------

         
library(ggplot2)        
library(gratia) 
library(marginaleffects) 
library(mvgam) 
library(parallel)
library(RColorBrewer)
library(stringr)
library(tibble)
library(tidyr)
library(tidyverse)    
library(wader)



everglades_counts_all <- tibble(max_counts(level = "all"))
everglades_counts_all <- everglades_counts_all |>
  filter(species %in% c("gbhe", "greg", "sneg", "whib", "wost", 'rosp')) |>  
  complete(year = full_seq(year, 1), species, fill = list(count = 0)) |>
  ungroup()




# setup plots  ------------------------------------------------------------
#can edit this later, but like how they look 


theme_set(theme_classic(base_size = 12, base_family = 'serif') +
            theme(axis.line.x.bottom = element_line(colour = "black",
                                                    size = 1),
                  axis.line.y.left = element_line(colour = "black",
                                                  size = 1)))

#display.brewer.all() 
#add different colors


options(ggplot2.discrete.colour = 'black',
        ggplot2.discrete.fill = brewer.pal(9, name = 'Paired'))


# add water data ----------------------------------------------------------


water <- load_datafile("eden_covariates.csv")
everglades_water_all <- filter(water, region == "all") 

count_env_data_all <- everglades_counts_all |>
  filter(year >= 1991) |> # No water data prior to 1991
  full_join(everglades_water_all, by = "year") |>
  drop_na(species) |>
  mutate(time = year - min(year) + 1, 
         series = factor(species) ) 




get_mvgam_priors(count ~ 1,
                 data = count_env_data_all,
                 family = gaussian()
)




# split data into training and testing sets -------------------------------


data_train <- filter(count_env_data_all, year < 2022 )
data_test <- filter(count_env_data_all, year >= 2022 )


# priors ------------------------------------------------------------------

def_priors <- get_mvgam_priors(
  # Observation formula containing species intercepts
  formula = count ~ series,
  
  # Process model that contains the hierarchical temporal smooths
  trend_formula = ~
    0 +
    
    # Shared smooth of year for all series
    s(year, 
      # there are 13 years in the training data
      k = 13, 
      bs = 'cr') +
    
    # Deviation smooths for each species
    s(year, 
      trend, 
      # there are 13 years in the training data
      k = 13, 
      bs = 'sz', 
      xt = list(bs = 'cr')),
  data = data_train,
  family = nb()
)

def_priors

sigma_prior <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)
intercept_prior <- prior(normal(0, 0.001), class = Intercept)
ar_sp_intercept_prior <- prior(std_normal(), class = b)

gam_var_priors <- c(sigma_prior, ar_sp_intercept_prior, intercept_prior)


# Models ------------------------------------------------------------------

#Null - mean
baseline_model_mean <- mvgam( formula = count ~ series, 
                    data = data_train,
                    newdata = data_test
)

null_AR_mvgam <- mvgam( formula = count ~ 1, 
                        treand_model = AR(), 
                          data = data_train,
                          newdata = data_test
)
null_VAR_mvgam <- mvgam( formula = count ~ 1, 
                         trend_model = VAR(), 
                          data = data_train,
                          newdata = data_test
)


#Trait
trait_mvgam <- mvgam(
  # Observation formula containing species-level intercepts
  formula = count ~ series,
  
  # Process model that contains the hierarchical temporal smooths
  trend_formula = ~
    0 + #think adding 0 makes it structured. 
    
    # Shared smooth of x for all series
    s(init_depth, 
      bs = 'cr') +
    
    # Deviation smooths for each series
    s(init_depth,
      trend,
      bs = 'sz',
      xt = list(bs = 'cr'))+
    
    # Shared smooth of x for all series
    s(dry_days, 
      bs = 'cr') +
    

    
    # Shared smooth of x for all series
    s(breed_season_depth, 
      bs = 'cr') +
    # 
    # # Deviation smooths for each series
    s(breed_season_depth,
      trend,
      bs = 'sz',
      xt = list(bs = 'cr'))+
    
    
    # Shared smooth of x for all series
    s(recession, 
      bs = 'cr') #+
  
  # # Deviation smooths for each series
  # s(recession,                                       #removing this makes R_hat larger
  #   trend,
  #   bs = 'sz',
  #   xt = list(bs = 'cr'))
  ,
  
  trend_model = VAR(), 
  
  trend_map =
    # trend_map forces species to track same /different latent signals
    data.frame(
      series = unique(data_train$series),
      trend = c(1, 1, 2, 1, 3, 4)
    ),

  
  # Updated prior distributions for the series-level 
  # intercepts using brms::prior()
  priors = ar_sp_intercept_prior,
  
  # Training and testing data in mvgam's long format
  data = data_train,
  newdata = data_test,
  
  # nb observation model
  family = nb(),
  
  
  control = list(max_treedepth = 10,   #10 
                 adapt_delta = 0.9),   #0.8
  
  # Each series shares the same nb shape parameter
  # If TRUE and the family has additional 
  # family-specific observation parameters (e.g., 
  # variance components, dispersion parameters), 
  # these will be shared across all outcome variables. 
  # Useful when multiple outcomes share properties. 
  share_obs_params = TRUE,
  
  # Non-centring the latent states tends to improve
  # performance in State Space models
  noncentred = TRUE,
  backend = 'cmdstanr',
  
  
  
  samples = 1500
)



trait_mvgam_variables <- mvgam(
  # Observation formula containing species-level intercepts
  formula = count ~ series,
  
  # Process model that contains the hierarchical temporal smooths
  trend_formula = ~
    0 + 
    s(init_depth,
      trend,
      bs = 'sz',
      xt = list(bs = 'cr'))+
    s(dry_days, 
      bs = 'cr') +
    s(breed_season_depth,
      trend,
      bs = 'sz',
      xt = list(bs = 'cr'))+
    dry_days 
  ,
  
  trend_model = VAR(), 
  
  trend_map =
    # trend_map forces species to track same /different latent signals
    data.frame(
      series = unique(data_train$series),
      trend = c(1, 1, 2, 1, 3, 2)
    ),

  
  
  # Updated prior distributions for the series-level 
  # intercepts using brms::prior()
  priors = ar_sp_intercept_prior,
  
  # Training and testing data in mvgam's long format
  data = data_train,
  newdata = data_test,
  
  # nb observation model
  family = nb(),
  
  
  control = list(max_treedepth = 10,   #10 
                 adapt_delta = 0.9),   #0.8
  
  # Each series shares the same nb shape parameter
  # If TRUE and the family has additional 
  # family-specific observation parameters (e.g., 
  # variance components, dispersion parameters), 
  # these will be shared across all outcome variables. 
  # Useful when multiple outcomes share properties. 
  share_obs_params = TRUE,
  
  # Non-centring the latent states tends to improve
  # performance in State Space models
  noncentred = TRUE,
  backend = 'cmdstanr',
  
  
  
  samples = 1500
)


var_mvgam <- mvgam(
  formula = count ~ 1,
  trend_formula = ~ s(breed_season_depth, trend, bs = "re") +
    s(dry_days, trend, bs = "re") +
    s(recession, trend, bs = "re"),
  trend_model = VAR(),
  family = nb(),
  data = data_train,
  newdata = data_test,
  chains = 4
)

ar_mvgam <- mvgam(
  formula = count ~ 1,
  trend_formula = ~ s(breed_season_depth, trend, bs = "re") +
    s(dry_days, trend, bs = "re") +
    s(recession, trend, bs = "re"),
  trend_model = AR(),
  family = nb(),
  data = data_train,
  newdata = data_test,
  chains = 4
)






trait_mvgam_variables2 <- mvgam(
  # Observation formula containing species-level intercepts
  formula = count ~ series,
  
  # Process model that contains the hierarchical temporal smooths
  trend_formula = ~
    0 + 
    s(init_depth,
      trend,
      bs = 'sz',
      xt = list(bs = 'cr'))+
    s(dry_days, 
      bs = 'cr') +
    s(breed_season_depth,
      trend,
      bs = 'sz',
      xt = list(bs = 'cr'))+
    dry_days 
  ,
  
  trend_model = AR(), 
  
  trend_map =
    # trend_map forces species to track same /different latent signals
    data.frame(
      series = unique(data_train$series),
      trend = c(1, 1, 2, 1, 2, 2)
    ),
  

  
  # Updated prior distributions for the series-level 
  # intercepts using brms::prior()
  priors = ar_sp_intercept_prior,
  
  # Training and testing data in mvgam's long format
  data = data_train,
  newdata = data_test,
  
  # nb observation model
  family = nb(),
  
  
  control = list(max_treedepth = 10,   #10 
                 adapt_delta = 0.9),   #0.8
  
  # Each series shares the same nb shape parameter
  # If TRUE and the family has additional 
  # family-specific observation parameters (e.g., 
  # variance components, dispersion parameters), 
  # these will be shared across all outcome variables. 
  # Useful when multiple outcomes share properties. 
  share_obs_params = TRUE,
  
  # Non-centring the latent states tends to improve
  # performance in State Space models
  noncentred = TRUE,
  backend = 'cmdstanr',
  
  
  
  samples = 1500
)


var_mvgam <- mvgam(
  formula = count ~ 1,
  trend_formula = ~ s(breed_season_depth, trend, bs = "re") +
    s(dry_days, trend, bs = "re") +
    s(recession, trend, bs = "re"),
  trend_model = VAR(),
  family = nb(),
  data = data_train,
  newdata = data_test,
  chains = 4
)

ar_mvgam <- mvgam(
  formula = count ~ 1,
  trend_formula = ~ s(breed_season_depth, trend, bs = "re") +
    s(dry_days, trend, bs = "re") +
    s(recession, trend, bs = "re"),
  trend_model = AR(),
  family = nb(),
  data = data_train,
  newdata = data_test,
  chains = 4
)


# Forecasting - plotting --------------------------------------------------


#for later

#Null
fc_mean_null <- forecast(
  baseline_model_mean)

fc_AR_null <- forecast(
  null_AR_mvgam)

fc_VAR_null <- forecast(
  null_VAR_mvgam)

# forecast_null1 <- plot(fc_mean_null, series = 1,
#                        main = "Null model")
#  forecast_null2 <- plot(fc_mean_null, series = 2,
#                        main = "Null model")
#  forecast_null3 <- plot(fc_mean_null, series = 3,
#                        main = "Null model")
#  forecast_null4 <- plot(fc_mean_null, series = 4,
#                        main = "Null model")
#  forecast_null5 <- plot(fc_mean_null, series = 5,
#                        main = "Null model")
#  forecast_null6 <- plot(fc_mean_null, series = 6,
#                        main = "Null model")



#Trait 
fc_trait <- forecast(
  trait_mvgam)

# forecast_trait <- plot(fc_trait, series = 1, 
#      main = "Trait model")
# forecast_trait <- plot(fc_trait, series = 2, 
#      main = "Trait model")
# forecast_trait <- plot(fc_trait, series = 3, 
#      main = "Trait model")
# forecast_trait <- plot(fc_trait, series = 4, 
#      main = "Trait model")
# forecast_trait <- plot(fc_trait, series = 5, 
#      main = "Trait model")
# forecast_trait <- plot(fc_trait, series = 6, 
#      main = "Trait model")

fc_tv <-forecast(
  trait_mvgam_variables)








# scores ------------------------------------------------------------------


fc_mean_null <- forecast(baseline_model_mean, 
                    score = 'crps')
null_mean_score <- score(fc_mean_null, 
                    score = 'crps')
null_mean_score <- mapply(cbind, null_mean_score, "model"= 'null_mean', SIMPLIFY=F)


fc_AR_null <- forecast(null_AR_mvgam, 
                    score = 'crps')
null_AR_score <- score(fc_AR_null, 
                    score = 'crps')
null_AR_score <- mapply(cbind, null_AR_score, "model"= 'null_ar', SIMPLIFY=F)


fc_VAR_null <- forecast(null_VAR_mvgam, 
                    score = 'crps')
null_VAR_score <- score(fc_VAR_null, 
                    score = 'crps')
null_VAR_score <- mapply(cbind, null_VAR_score, "model"= 'null_var', SIMPLIFY=F)



fc_trait <- forecast(trait_mvgam, 
                    score = 'crps')
trait_score <- score(fc_trait, 
                    score = 'crps')
trait_score <- mapply(cbind, trait_score, "model"= 'trait', SIMPLIFY=F)



fc_tv <- forecast(trait_mvgam_variables, 
                     score = 'crps')
tv_score <- score(fc_tv, 
                     score = 'crps')
tv_score <- mapply(cbind, tv_score, "model"= 'tv', SIMPLIFY=F)



fc_tv2 <- forecast(trait_mvgam_variables2, 
                  score = 'crps')
tv2_score <- score(fc_tv2, 
                  score = 'crps')
tv2_score <- mapply(cbind, tv2_score, "model"= 'tv2', SIMPLIFY=F)


fc_var <- forecast(var_mvgam, 
                     score = 'crps')
var_score <- score(fc_var, 
                     score = 'crps')
var_score <- mapply(cbind, var_score, "model"= 'var', SIMPLIFY=F)



fc_ar <- forecast(ar_mvgam, 
                     score = 'crps')
ar_score <- score(fc_ar, 
                     score = 'crps')
ar_score <- mapply(cbind, ar_score, "model"= 'ar', SIMPLIFY=F)




scores_species_fig <- rbind(
  within(Map(cbind, null_mean_score, group = names(null_mean_score)), 
         rm(all_series)) |> 
    bind_rows(), 
  
  # within(Map(cbind, null_AR_score, group = names(null_AR_score)), 
  #        rm(all_series)) |> 
  #   bind_rows(), 
  
  # within(Map(cbind, null_VAR_score, group = names(null_VAR_score)), 
  #        rm(all_series)) |> 
  #   bind_rows(), 
  
  within(Map(cbind, trait_score, group = names(trait_score)),
         rm(all_series)) |>
    bind_rows(),
   
  # within(Map(cbind, var_score, group = names(var_score)), 
  #        rm(all_series)) |> 
  #   bind_rows(),
  
  within(Map(cbind, tv_score, group = names(tv_score)), 
         rm(all_series)) |> 
    bind_rows(),
  
  
  within(Map(cbind, tv2_score, group = names(tv2_score)), 
         rm(all_series)) |> 
    bind_rows()
  
  # within(Map(cbind, ar_score, group = names(ar_score)), 
  #        rm(all_series))
  # |> 
  #   bind_rows() 
)|> 
  arrange(group, eval_horizon, score) |> 
  ggplot(aes(x = eval_horizon, y = score, fill= model)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~group, scales = "free")



model_all_scores <- rbind(
  null_mean_score$all_series,
 trait_score$all_series,
 # null_AR_score$all_series,
 # null_VAR_score$all_series,
 # var_score$all_series, 
 tv_score$all_series,
 tv2_score$all_series
 ) |> 
  arrange(eval_horizon, score) |> 
  mutate(model = factor(model))



scores_all_fig <- model_all_scores |> 
  ggplot(aes(x = eval_horizon, y = score, fill= model)) +
  geom_bar(stat="identity", position="dodge")

scores_species_fig
scores_all_fig

#trait: VAR
# gbhe=greg=sneg
# rosp  
# whib 
# wost

#tv: VAR
# gbhe=greg=sneg 
# rosp=wost
# whib

#tv2:AR
# gbhe=greg=sneg 
# rosp=whib=wost


 fc_mean_null |> plot()
 fc_tv |> plot()
 fc_tv2 |> plot()
 fc_trait |> plot()

# scores_fig

ggsave('results/scores/scores_all_fig.png', scores_all_fig, 
       width = 12, height = 6, dpi = 300)


ggsave('results/scores/scores_species_fig.png', scores_species_fig, 
       width = 12, height = 6, dpi = 300)




#need trend model in my mvgams - for rolling forecast... 






# option to display scores ------------------------------------------------

# crps_mod1 <- score(fc_mod1, score = "crps")
# crps_mod2 <- score(fc_mod2, score = "crps")
# 
# diff_scores <- crps_mod2$series_1$score -
#   crps_mod1$series_1$score
# plot(diff_scores,
#      pch = 16, cex = 1.25, col = "darkred",
#      ylim = c(
#        -1 * max(abs(diff_scores), na.rm = TRUE),
#        max(abs(diff_scores), na.rm = TRUE)
#      ),
#      bty = "l",
#      xlab = "Forecast horizon",
#      ylab = expression(CRPS[AR1] ~ -~ CRPS[spline])
# )
# abline(h = 0, lty = "dashed", lwd = 2)
# ar1_better <- length(which(diff_scores < 0))
# title(main = paste0(
#   "AR(1) better in ",
#   ar1_better,
#   " of ",
#   length(diff_scores),
#   " evaluations",
#   "\nMean difference = ",
#   round(mean(diff_scores, na.rm = TRUE), 2)
# ))

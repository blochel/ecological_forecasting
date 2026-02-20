
# DGAM models - Everglades Wading Birds -----------------------------------



library(mvgam)          
library(tidyverse)      
library(ggplot2)         
library(marginaleffects) 
library(gratia)         
library(parallel)
library(stringr)
library(tibble)
library(tidyr)
library(wader)



everglades_counts_all <- tibble(max_counts(level = "all"))
everglades_counts_all <- everglades_counts_all |>
  filter(species %in% c("gbhe", "greg", "sneg", "whib", "wost", 'rosp')) |>  
  complete(year = full_seq(year, 1), species, fill = list(count = 0)) |>
  ungroup()





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

#Null 
null_mvgam <- mvgam(count ~
                      1, 
                    trend_model = AR(),
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
    
    # Shared smooth of x for all series
    s(dry_days, 
      bs = 'cr') +
    
    # Deviation smooths for each series
    s(dry_days,
      trend,
      bs = 'sz',
      xt = list(bs = 'cr'))+
    
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




# Forecasting - plotting --------------------------------------------------


#for later

#Null
fc_null <- forecast(
  null_mvgam)

# forecast_null <- plot(fc_null, series = 1, 
#                       main = "Null model")
# forecast_null <- plot(fc_null, series = 2, 
#                       main = "Null model")
# forecast_null <- plot(fc_null, series = 3, 
#                       main = "Null model")
# forecast_null <- plot(fc_null, series = 4, 
#                       main = "Null model")
# forecast_null <- plot(fc_null, series = 5, 
#                       main = "Null model")
# forecast_null <- plot(fc_null, series = 6, 
#                       main = "Null model")



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












# scores ------------------------------------------------------------------


fc_null <- forecast(null_mvgam, 
                    score = 'crps')
null_score <- score(fc_null, 
                    score = 'crps')
null_score <- mapply(cbind, null_score, "model"= 'null', SIMPLIFY=F)




fc_trait <- forecast(trait_mvgam, 
                    score = 'crps')
trait_score <- score(fc_trait, 
                    score = 'crps')
trait_score <- mapply(cbind, trait_score, "model"= 'trait', SIMPLIFY=F)





model_scores <- rbind(
  trait_score$all_series,
  null_score$all_series) |> 
  arrange(eval_horizon, score) |> 
  mutate(model = factor(model))



scores_fig <- model_scores |> 
  ggplot(aes(x = eval_horizon, y = score, fill= model)) +
  geom_bar(stat="identity", position="dodge")

# scores_fig

ggsave('results/scores/scores_fig.png', scores_fig, width = 12, height = 6, dpi = 300)

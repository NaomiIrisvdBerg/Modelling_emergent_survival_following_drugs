##Fitting monoculture drug sensitivities 

##use objects defined in other script ('Bug_drug_community_effects_Sept24.R', first ~2100 lines:
#The object below was made from the object 'combined_monoculture_aucs' but has a control row (AUC == 1) for all drugs
combined_monoculture_aucs1
#The object below is made from the object 'combined_monoculture_aucs1' and contains drug dynamic info per drug as well (as deduced from metabolomics data)
combined_monoculture_aucs2

str(combined_monoculture_aucs2)

#load required packages
library(dplyr)
library(ggplot2)
library(vegan)
library(minpack.lm)
library(dplyr)

##Ensure all unique combination of species and drugs have their control AUC (1) in the dataset as well
species_drug_combinations <- combined_monoculture_aucs2 %>%
  dplyr::group_by(Species, species_name, drug, NT_code, Prop_left_24h_background, cc_K) %>%
  dplyr::filter(all(conc > 0)) %>%
  dplyr::distinct(Species, NT_code, drug)

new_rows <- species_drug_combinations %>%
  dplyr::mutate(conc = 0, 
                concentration = 0, 
                effective_concentration = 0, 
                AUC = 1,
                control = TRUE, 
                abs_change_AUC = 0,
                AUC_cc_K = NA)  ##not important in fitting since we will use AUCs directly


combined_monoculture_aucs3 <- bind_rows(combined_monoculture_aucs2, new_rows)


combined_monoculture_aucs3_Acarbose <- subset(combined_monoculture_aucs3, drug == "Acarbose")
combined_monoculture_aucs3_Amlodipine <- subset(combined_monoculture_aucs3, drug == "Amlodipine")
combined_monoculture_aucs3_Aprepitant <- subset(combined_monoculture_aucs3, drug == "Aprepitant")
combined_monoculture_aucs3_Aripiprazole <- subset(combined_monoculture_aucs3, drug == "Aripiprazole")
combined_monoculture_aucs3_Azithromycin <- subset(combined_monoculture_aucs3, drug == "Azithromycin")
combined_monoculture_aucs3_beta_estradiol_17_valerate <- subset(combined_monoculture_aucs3, drug == "beta-estradiol 17-valerate")
combined_monoculture_aucs3_Chlorpromazine <- subset(combined_monoculture_aucs3, drug == "Chlorpromazine")
combined_monoculture_aucs3_Ciprofloxacin <- subset(combined_monoculture_aucs3, drug == "Ciprofloxacin")
combined_monoculture_aucs3_Clomiphene <- subset(combined_monoculture_aucs3, drug == "Clomiphene")
combined_monoculture_aucs3_Diacerein <- subset(combined_monoculture_aucs3, drug == "Diacerein")
combined_monoculture_aucs3_Dienestrol <- subset(combined_monoculture_aucs3, drug == "Dienestrol")
combined_monoculture_aucs3_Doxycycline <- subset(combined_monoculture_aucs3, drug == "Doxycycline")
combined_monoculture_aucs3_Ebselen <- subset(combined_monoculture_aucs3, drug == "Ebselen")
combined_monoculture_aucs3_Entacapone <- subset(combined_monoculture_aucs3, drug == "Entacapone")
combined_monoculture_aucs3_Felodipine <- subset(combined_monoculture_aucs3, drug == "Felodipine")
combined_monoculture_aucs3_Ketoconazole <- subset(combined_monoculture_aucs3, drug == "Ketoconazole")
combined_monoculture_aucs3_Lansoprazole <- subset(combined_monoculture_aucs3, drug == "Lansoprazole")
combined_monoculture_aucs3_Loratadine <- subset(combined_monoculture_aucs3, drug == "Loratadine")
combined_monoculture_aucs3_Loxapine <- subset(combined_monoculture_aucs3, drug == "Loxapine")
combined_monoculture_aucs3_Mecillinam <- subset(combined_monoculture_aucs3, drug == "Mecillinam")
combined_monoculture_aucs3_Mefloquine <- subset(combined_monoculture_aucs3, drug == "Mefloquine")
combined_monoculture_aucs3_Methotrexate <- subset(combined_monoculture_aucs3, drug == "Methotrexate")
combined_monoculture_aucs3_Niclosamide <- subset(combined_monoculture_aucs3, drug == "Niclosamide")
combined_monoculture_aucs3_Nifurtimox <- subset(combined_monoculture_aucs3, drug == "Nifurtimox")
combined_monoculture_aucs3_Olanzapine <- subset(combined_monoculture_aucs3, drug == "Olanzapine")
combined_monoculture_aucs3_Omeprazole <- subset(combined_monoculture_aucs3, drug == "Omeprazole")
combined_monoculture_aucs3_Oxolinic_acid <- subset(combined_monoculture_aucs3, drug == "Oxolinic acid")
combined_monoculture_aucs3_Sertindole <- subset(combined_monoculture_aucs3, drug == "Sertindole")
combined_monoculture_aucs3_Simvastatin <- subset(combined_monoculture_aucs3, drug == "Simvastatin")
combined_monoculture_aucs3_Tamoxifen <- subset(combined_monoculture_aucs3, drug == "Tamoxifen")

#where x = "effective_concentration"
#where y = "AUC"
#where y_at_max_x = "AUC" at the maximum measured "effective_concentration"

#where e and b are fitted parameters, aiming to have the best fit (minimal RSS) describing the change in "AUC" as a function of "effective_concentration"




####
##Acarbose
results_df_fits_Acarbose <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Acarbose$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Acarbose %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  #check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      d <- d
      d1 <- d1
      
      
      #Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Acarbose <- rbind(results_df_fits_Acarbose, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Acarbose <- rbind(results_df_fits_Acarbose, fit_results_df)
    }
  }
}

#View(results_df_fits_Acarbose)


##Amlodipine
results_df_fits_Amlodipine <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Amlodipine$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Amlodipine %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      d <- d
      d1 <- d1
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Amlodipine <- rbind(results_df_fits_Amlodipine, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Amlodipine <- rbind(results_df_fits_Amlodipine, fit_results_df)
    }
  }
}

#View(results_df_fits_Amlodipine)





##Aprepitant
results_df_fits_Aprepitant <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Aprepitant$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Aprepitant %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      d <- d
      d1 <- d1
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Aprepitant <- rbind(results_df_fits_Aprepitant, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Aprepitant <- rbind(results_df_fits_Aprepitant, fit_results_df)
    }
  }
}

#View(results_df_fits_Aprepitant)





##Aripiprazole
results_df_fits_Aripiprazole <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Aripiprazole$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Aripiprazole %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      d <- d
      d1 <- d1
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Aripiprazole <- rbind(results_df_fits_Aripiprazole, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Aripiprazole <- rbind(results_df_fits_Aripiprazole, fit_results_df)
    }
  }
}

#View(results_df_fits_Aripiprazole)






##Azithromycin
results_df_fits_Azithromycin <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Azithromycin$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Azithromycin %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      d <- d
      d1 <- d1
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Azithromycin <- rbind(results_df_fits_Azithromycin, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Azithromycin <- rbind(results_df_fits_Azithromycin, fit_results_df)
    }
  }
}

#View(results_df_fits_Azithromycin)





##beta_estradiol_17_valerate
results_df_fits_beta_estradiol_17_valerate <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_beta_estradiol_17_valerate$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_beta_estradiol_17_valerate %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  #Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      d <- d
      d1 <- d1
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_beta_estradiol_17_valerate <- rbind(results_df_fits_beta_estradiol_17_valerate, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_beta_estradiol_17_valerate <- rbind(results_df_fits_beta_estradiol_17_valerate, fit_results_df)
    }
  }
}

#View(results_df_fits_beta_estradiol_17_valerate)





##Chlorpromazine
results_df_fits_Chlorpromazine <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Chlorpromazine$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Chlorpromazine %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Chlorpromazine <- rbind(results_df_fits_Chlorpromazine, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Chlorpromazine <- rbind(results_df_fits_Chlorpromazine, fit_results_df)
    }
  }
}

#View(results_df_fits_Chlorpromazine)





##Ciprofloxacin
results_df_fits_Ciprofloxacin <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Ciprofloxacin$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Ciprofloxacin %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      d <- d
      d1 <- d1
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Ciprofloxacin <- rbind(results_df_fits_Ciprofloxacin, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Ciprofloxacin <- rbind(results_df_fits_Ciprofloxacin, fit_results_df)
    }
  }
}

#View(results_df_fits_Ciprofloxacin)




##Clomiphene
results_df_fits_Clomiphene <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Clomiphene$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Clomiphene %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Clomiphene <- rbind(results_df_fits_Clomiphene, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Clomiphene <- rbind(results_df_fits_Clomiphene, fit_results_df)
    }
  }
}

#View(results_df_fits_Clomiphene)






##Diacerein
results_df_fits_Diacerein <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Diacerein$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Diacerein %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Diacerein <- rbind(results_df_fits_Diacerein, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Diacerein <- rbind(results_df_fits_Diacerein, fit_results_df)
    }
  }
}

#View(results_df_fits_Diacerein)




##Dienestrol
results_df_fits_Dienestrol <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Dienestrol$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Dienestrol %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Dienestrol <- rbind(results_df_fits_Dienestrol, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Dienestrol <- rbind(results_df_fits_Dienestrol, fit_results_df)
    }
  }
}

#View(results_df_fits_Dienestrol)




##Doxycycline
results_df_fits_Doxycycline <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Doxycycline$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Doxycycline %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Doxycycline <- rbind(results_df_fits_Doxycycline, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Doxycycline <- rbind(results_df_fits_Doxycycline, fit_results_df)
    }
  }
}

#View(results_df_fits_Doxycycline)





##Ebselen
results_df_fits_Ebselen <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Ebselen$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Ebselen %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Ebselen <- rbind(results_df_fits_Ebselen, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Ebselen <- rbind(results_df_fits_Ebselen, fit_results_df)
    }
  }
}

#View(results_df_fits_Ebselen)




##Entacapone
results_df_fits_Entacapone <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Entacapone$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Entacapone %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Entacapone <- rbind(results_df_fits_Entacapone, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Entacapone <- rbind(results_df_fits_Entacapone, fit_results_df)
    }
  }
}

#View(results_df_fits_Entacapone)





##Felodipine
results_df_fits_Felodipine <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Felodipine$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Felodipine %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Felodipine <- rbind(results_df_fits_Felodipine, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Felodipine <- rbind(results_df_fits_Felodipine, fit_results_df)
    }
  }
}

#View(results_df_fits_Felodipine)





##Ketoconazole
results_df_fits_Ketoconazole <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Ketoconazole$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Ketoconazole %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Ketoconazole <- rbind(results_df_fits_Ketoconazole, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Ketoconazole <- rbind(results_df_fits_Ketoconazole, fit_results_df)
    }
  }
}

#View(results_df_fits_Ketoconazole)





##Lansoprazole
results_df_fits_Lansoprazole <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Lansoprazole$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Lansoprazole %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Lansoprazole <- rbind(results_df_fits_Lansoprazole, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Lansoprazole <- rbind(results_df_fits_Lansoprazole, fit_results_df)
    }
  }
}

#View(results_df_fits_Lansoprazole)



##Loratadine
results_df_fits_Loratadine <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Loratadine$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Loratadine %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Loratadine <- rbind(results_df_fits_Loratadine, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Loratadine <- rbind(results_df_fits_Loratadine, fit_results_df)
    }
  }
}

#View(results_df_fits_Loratadine)



##Loxapine
results_df_fits_Loxapine <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Loxapine$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Loxapine %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Loxapine <- rbind(results_df_fits_Loxapine, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Loxapine <- rbind(results_df_fits_Loxapine, fit_results_df)
    }
  }
}

#View(results_df_fits_Loxapine)





##Mecillinam
results_df_fits_Mecillinam <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Mecillinam$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Mecillinam %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Mecillinam <- rbind(results_df_fits_Mecillinam, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA, 
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Mecillinam <- rbind(results_df_fits_Mecillinam, fit_results_df)
    }
  }
}

#View(results_df_fits_Mecillinam)





##Mefloquine
results_df_fits_Mefloquine <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Mefloquine$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Mefloquine %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Mefloquine <- rbind(results_df_fits_Mefloquine, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Mefloquine <- rbind(results_df_fits_Mefloquine, fit_results_df)
    }
  }
}

#View(results_df_fits_Mefloquine)




##Methotrexate
results_df_fits_Methotrexate <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Methotrexate$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Methotrexate %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Methotrexate <- rbind(results_df_fits_Methotrexate, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Methotrexate <- rbind(results_df_fits_Methotrexate, fit_results_df)
    }
  }
}

#View(results_df_fits_Methotrexate)




##Niclosamide
results_df_fits_Niclosamide <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Niclosamide$Species)) {
  
  ##Subset data for the current species, and subset for concentrations tested in co-culture, just like the other drugs 
  #(Niclosamide is an exception since it has more data available, but for consistency we approach it like we do the other drugs)
  species_data <- subset(combined_monoculture_aucs3_Niclosamide, concentration %in% c(0, 20, 80, 160))  %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.1, b = 0.1), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.1, e = 0.1), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.1, e = 0.1), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Niclosamide <- rbind(results_df_fits_Niclosamide, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Niclosamide <- rbind(results_df_fits_Niclosamide, fit_results_df)
    }
  }
}

#View(results_df_fits_Niclosamide)



##Nifurtimox
results_df_fits_Nifurtimox <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Nifurtimox$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Nifurtimox %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Nifurtimox <- rbind(results_df_fits_Nifurtimox, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Nifurtimox <- rbind(results_df_fits_Nifurtimox, fit_results_df)
    }
  }
}

#View(results_df_fits_Nifurtimox)




##Olanzapine
results_df_fits_Olanzapine <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Olanzapine$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Olanzapine %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Olanzapine <- rbind(results_df_fits_Olanzapine, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Olanzapine <- rbind(results_df_fits_Olanzapine, fit_results_df)
    }
  }
}

#View(results_df_fits_Olanzapine)




##Omeprazole
results_df_fits_Omeprazole <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Omeprazole$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Omeprazole %>%
    dplyr::filter(Species == species)
  
  # Define x and y
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Omeprazole <- rbind(results_df_fits_Omeprazole, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Omeprazole <- rbind(results_df_fits_Omeprazole, fit_results_df)
    }
  }
}

#View(results_df_fits_Omeprazole)




##Oxolinic_acid
results_df_fits_Oxolinic_acid <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Oxolinic_acid$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Oxolinic_acid %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 2.1577658, e = 0.4514893), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 2.116853, e = 0.5336618), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Oxolinic_acid <- rbind(results_df_fits_Oxolinic_acid, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Oxolinic_acid <- rbind(results_df_fits_Oxolinic_acid, fit_results_df)
    }
  }
}

#View(results_df_fits_Oxolinic_acid)




##Sertindole
results_df_fits_Sertindole <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Sertindole$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Sertindole %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Sertindole <- rbind(results_df_fits_Sertindole, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Sertindole <- rbind(results_df_fits_Sertindole, fit_results_df)
    }
  }
}

#View(results_df_fits_Sertindole)




##Simvastatin
results_df_fits_Simvastatin <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Simvastatin$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Simvastatin %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) *exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA)
    
  )
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Simvastatin <- rbind(results_df_fits_Simvastatin, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Simvastatin <- rbind(results_df_fits_Simvastatin, fit_results_df)
    }
  }
}

#View(results_df_fits_Simvastatin)






##Tamoxifen
results_df_fits_Tamoxifen <- data.frame()

##Loop through each unique species
for (species in unique(combined_monoculture_aucs3_Tamoxifen$Species)) {
  
  ##Subset data for the current species
  species_data <- combined_monoculture_aucs3_Tamoxifen %>%
    dplyr::filter(Species == species)
  
  ##Define x and y and other predetermined params
  x <- species_data$effective_concentration
  max_x <- max(x)
  y <- species_data$AUC
  #d <- y[which.max(x)]
  d1 <- max(y, na.rm = TRUE)
  
  sorted_x <- sort(x, decreasing = TRUE)
  second_max_x <- sorted_x[2]
  second_max_y <- y[x == second_max_x]
  
  d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])
  
  ##Adjust x values to avoid fitting issues at zero
  x <- ifelse(x == 0, x + 0.0001, x)
  
  ##Check if there are enough data points to fit
  if(length(x) < 3) next
  
  ##Define all models and fit them
  models <- list(
    Linear_func = tryCatch(nlsLM(y ~ 1 + e * x, start = list(e = 0.01), 
                                 control = nls.lm.control(maxiter = 1024)),
                           error = function(e) NA),
    Asymptotic_func = tryCatch(nlsLM(y ~ (d - (d - 1) * exp(-e * x)),
                                     start = list(e = 1), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Log_log_func = tryCatch(nlsLM(y ~ 1 + ((d - 1) / (1 + exp(-e * (log10(x) - log10(b))))),
                                  start = list(e = 0.01, b = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-b * (x - e))),
                                   start = list(b = 0.01, e = 0.01), 
                                   control = nls.lm.control(maxiter = 1024)),
                             error = function(e) NA),
    Modified_Gompertz_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(b*(x-e)))),
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Inverse_sigmoidal_func = tryCatch(nlsLM(y ~ ((1 - d) / (1 + exp(b * x - e))) + d,
                                            start = list(b = 0.01, e = 0.01), 
                                            control = nls.lm.control(maxiter = 1024)),
                                      error = function(e) NA),
    Weibull_I_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * exp(-exp(-e * (log10(x) - log10(b)))),
                                    start = list(e = 0.01, b = 0.01), 
                                    control = nls.lm.control(maxiter = 1024)),
                              error = function(e) NA),
    Weibull_II_func = tryCatch(nlsLM(y ~ 1 + (d - 1) * (1 - exp(-exp(e * (log10(x) - log10(b))))),
                                     start = list(e = 0.01, b = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Lorentz_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) / (1 + e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Lorentz_func_II = tryCatch(nlsLM(y ~ d + (d1 - d) / (1 + e * (x - b)^2),
                                     start = list(b = 0.01, e = 0.01), 
                                     control = nls.lm.control(maxiter = 1024)),
                               error = function(e) NA),
    Bragg_func = tryCatch(nlsLM(y ~ 1 + (d1 - 1) * exp(-e*(x-b)^2),
                                start = list(b = 0.01, e = 0.01), 
                                control = nls.lm.control(maxiter = 1024)),
                          error = function(e) NA),
    Mod_exp_func = tryCatch(nlsLM(y ~ d1 * exp(-e * x) + b * x,
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA),
    Gauss_peak_func = tryCatch(nlsLM(y ~ d + (d1 - d) * exp(-e * (x - b)^2),
                                  start = list(b = 0.01, e = 0.01), 
                                  control = nls.lm.control(maxiter = 1024)),
                            error = function(e) NA)
    
  )
  
  
  ##Fit each model and collect the parameters and RSS
  for (model_name in names(models)) {
    fit <- models[[model_name]]
    if (inherits(fit, "nls")) {  ##Check if the fit was successful
      params <- coef(fit)
      rss <- sum(resid(fit)^2)
      d <- d
      d1 <- d1
      
      ##Create a dataframe for the current fit results
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = rss,
        d = as.numeric(d),
        d1 = as.numeric(d1),
        Fit_param_name = names(params),
        Fit_param_value = as.numeric(params)
      )
      
      ##Append the current fit results to the main results dataframe
      results_df_fits_Tamoxifen <- rbind(results_df_fits_Tamoxifen, fit_results_df)
    } else {
      ##If fit failed, add an entry with NA values
      fit_results_df <- data.frame(
        Model_type = model_name,
        Species = species,
        Sum_of_squares = NA,
        d = NA,
        d1 = NA,
        Fit_param_name = NA,
        Fit_param_value = NA
      )
      
      ##Append the failed fit results to the main results dataframe
      results_df_fits_Tamoxifen <- rbind(results_df_fits_Tamoxifen, fit_results_df)
    }
  }
}

#View(results_df_fits_Tamoxifen)


results_df_fits_Acarbose$drug <- "Acarbose"
results_df_fits_Amlodipine$drug <- "Amlodipine"
results_df_fits_Aprepitant$drug <- "Aprepitant"
results_df_fits_Aripiprazole$drug <- "Aripiprazole"
results_df_fits_Azithromycin$drug <- "Azithromycin"
results_df_fits_beta_estradiol_17_valerate$drug <- "beta-estradiol 17-valerate"
results_df_fits_Chlorpromazine$drug <- "Chlorpromazine"
results_df_fits_Ciprofloxacin$drug <- "Ciprofloxacin"
results_df_fits_Clomiphene$drug <- "Clomiphene"
results_df_fits_Diacerein$drug <- "Diacerein"
results_df_fits_Dienestrol$drug <- "Dienestrol"
results_df_fits_Doxycycline$drug <- "Doxycycline"
results_df_fits_Ebselen$drug <- "Ebselen"
results_df_fits_Entacapone$drug <- "Entacapone"
results_df_fits_Felodipine$drug <- "Felodipine"
results_df_fits_Ketoconazole$drug <- "Ketoconazole"
results_df_fits_Lansoprazole$drug <- "Lansoprazole"
results_df_fits_Loratadine$drug <- "Loratadine"
results_df_fits_Loxapine$drug <- "Loxapine"
results_df_fits_Mecillinam$drug <- "Mecillinam"
results_df_fits_Mefloquine$drug <- "Mefloquine"
results_df_fits_Methotrexate$drug <- "Methotrexate"
results_df_fits_Niclosamide$drug <- "Niclosamide"
results_df_fits_Nifurtimox$drug <- "Nifurtimox"
results_df_fits_Olanzapine$drug <- "Olanzapine"
results_df_fits_Omeprazole$drug <- "Omeprazole"
results_df_fits_Oxolinic_acid$drug <- "Oxolinic acid"
results_df_fits_Sertindole$drug <- "Sertindole"
results_df_fits_Simvastatin$drug <- "Simvastatin"
results_df_fits_Tamoxifen$drug <- "Tamoxifen"



list_of_fits <- list(
  results_df_fits_Acarbose,
  results_df_fits_Amlodipine,
  results_df_fits_Aprepitant,
  results_df_fits_Aripiprazole,
  results_df_fits_Azithromycin,
  results_df_fits_beta_estradiol_17_valerate,
  results_df_fits_Chlorpromazine,
  results_df_fits_Ciprofloxacin,
  results_df_fits_Clomiphene,
  results_df_fits_Diacerein,
  results_df_fits_Dienestrol,
  results_df_fits_Doxycycline,
  results_df_fits_Ebselen,
  results_df_fits_Entacapone,
  results_df_fits_Felodipine,
  results_df_fits_Ketoconazole,
  results_df_fits_Lansoprazole,
  results_df_fits_Loratadine,
  results_df_fits_Loxapine,
  results_df_fits_Mecillinam,
  results_df_fits_Mefloquine,
  results_df_fits_Methotrexate,
  results_df_fits_Niclosamide,
  results_df_fits_Nifurtimox,
  results_df_fits_Olanzapine,
  results_df_fits_Omeprazole,
  results_df_fits_Oxolinic_acid,
  results_df_fits_Sertindole,
  results_df_fits_Simvastatin,
  results_df_fits_Tamoxifen
)

combined_fits_all_df <- do.call(rbind, list_of_fits)
#write.csv(combined_fits_all_df, file = "~/Downloads/combined_fits_all_df.csv", row.names = FALSE)


best_fitting_models_df <- combined_fits_all_df %>%
  dplyr::group_by(Species, drug) %>%
  dplyr::filter(Sum_of_squares == min(Sum_of_squares, na.rm = TRUE)) %>%
  dplyr::ungroup()


#best_fitting_models_df_Acarbose <- subset(best_fitting_models_df, drug == "Acarbose")



##take out '.' from species formatting to match with model formatting of species names/state variables:

best_fitting_models_df1 <- best_fitting_models_df
best_fitting_models_df1$Species <- gsub("\\.", "", best_fitting_models_df1$Species)
str(best_fitting_models_df1$Species)


#Here we have the numeric output:
# Function to generate y(Drug) with numeric values embedded
get_y_function <- function(species, drug, best_fitting_models_df1) {
  # Subset the dataframe for the specific species and drug
  model_info <- best_fitting_models_df1 %>%
    filter(Species == species, drug == drug)
  
  # Extract the model type
  model_type <- unique(model_info$Model_type)
  
  # Extract the parameters
  d <- unique(model_info$d)
  d1 <- unique(model_info$d1)
  params <- setNames(model_info$Fit_param_value, model_info$Fit_param_name)
  
  # Create the function string based on the model type
  function_string <- switch(model_type,
                            "Linear_func" = sprintf("function(Drug) 1 + %f * Drug", params["e"]),
                            "Asymptotic_func" = sprintf("function(Drug) %f - (%f - 1) * exp(-%f * Drug)", d, d, params["e"]),
                            "Log_log_func" = sprintf("function(Drug) 1 + ((%f - 1) / (1 + exp(-%f * (log10(Drug) - log10(%f)))))", d, params["e"], params["b"]),
                            "Gompertz_func" = sprintf("function(Drug) 1 + (%f - 1) * exp(-exp(-%f * (Drug - %f)))", d, params["b"], params["e"]),
                            "Modified_Gompertz_func" = sprintf("function(Drug) 1 + (%f - 1) * (1 - exp(-exp(%f*(Drug-%f))))", d, params["b"], params["e"]),
                            "Inverse_sigmoidal_func" = sprintf("function(Drug) ((1 - %f) / (1 + exp(%f * Drug - %f))) + %f", d, params["b"], params["e"], d),
                            "Weibull_I_func" = sprintf("function(Drug) 1 + (%f - 1) * exp(-exp(-%f * (log10(Drug) - log10(%f))))", d, params["e"], params["b"]),
                            "Weibull_II_func" = sprintf("function(Drug) 1 + (%f - 1) * (1 - exp(-exp(%f * (log10(Drug) - log10(%f)))))", d, params["e"], params["b"]),
                            "Lorentz_func" = sprintf("function(Drug) 1 + (%f - 1) / (1 + %f * (Drug - %f)^2)", d1, params["e"], params["b"]),
                            "Bragg_func" = sprintf("function(Drug) 1 + (%f - 1) * exp(-%f*(Drug-%f)^2)", d1, params["e"], params["b"]),
                            "Mod_exp_func" = sprintf("function(Drug) %f * exp(-%f * Drug) + %f * Drug", d1, params["e"], params["b"]),
                            "Gauss_peak_func" = sprintf("function(Drug) %f + (%f - %f) * exp(-%f * (Drug - %f)^2)", d, d1, d, params["e"], params["b"]),
                            "Lorentz_func_II" = sprintf("function(Drug) %f + (%f - %f) / (1 + %f * (Drug - %f)^2)", d, d1, d, params["e"], params["b"]),
                            stop("Model type not recognised")
  )
  
  
  #Parse and evaluate the function string to create the actual function
  y_function_eval <- eval(parse(text = function_string))
  
  return(y_function_eval)
}

#example usage
#species_elected <- "B._vulgatus"
#drug_elected <- "Acarbose"
#y_function <- get_y_function(species_elected, drug_elected, best_fitting_models_df_Acarbose)


best_fitting_models_df_Acarbose <- subset(best_fitting_models_df1, drug == "Acarbose")

y_function_B_vulgatus_Acarbose <- get_y_function("B_vulgatus", "Acarbose", best_fitting_models_df_Acarbose)

print(y_function_B_vulgatus_Acarbose)







####VISUALISATION OF FITS

##Plot those of interest
#Species formatting example: 'E._coli'
#Drugs capatalised except for beta_estradiol_17_valerate
#Drugs are a single word (no space), except for: beta_estradiol_17_valerate & Oxolinic_acid
species_elected <- 'E._coli'
drug_elected <- 'Chlorpromazine'

dataset_name <- paste0("combined_monoculture_aucs3_", drug_elected)

data_elected <- get(dataset_name)

#If Niclosamide: select data from minimally tested conc. (in co-cult) onwards (20 uM) to not bias rss to lower concentrations & make fitting more comparable to other drugs
#data_elected <- subset(data_elected,concentration %in% c(0, 20, 40, 80, 160))

concentrations <- as.numeric(data_elected$effective_concentration[data_elected$Species == species_elected])
x <- ifelse(concentrations == 0, concentrations + 0.0001, concentrations)

# Extract 'AUC' values for species selected
y <- data_elected$AUC[data_elected$Species == species_elected]

max_x <- max(x)       
#d <- max(y[x == max_x])
d1 <- max(y, na.rm = TRUE)

sorted_x <- sort(x, decreasing = TRUE)
second_max_x <- sorted_x[2]
second_max_y <- y[x == second_max_x]

d <- ifelse(is.na(y[x == max_x]), y[x == second_max_x], y[x == max_x])

plot(x, y, xlab= paste0("Effective concentration ", drug_elected), ylab='AUC', main = species_elected, pch=16, cex = 1.5, col = rgb(red = 0.5, green = 0, blue = 1, alpha = 0.3))


##Linear fit
Linear_func <- function(x, e){1+e*x}
bfit_linear <- nlsLM(y ~ Linear_func(x, e), start = list(e = 0))
bfit_linear
linear_coefs <- coef(bfit_linear)
linear_rss <- sum(resid(bfit_linear)^2)
linear_rss
curve(Linear_func(x, linear_coefs['e']), add = TRUE, col = "grey", lwd = 2)


##Asymptotic fit
asymptotic_func <- function(x, d, e) {
  (d - (d - 1) * exp(-e * x))
}

bfit_asymptotic <- nlsLM(y ~ asymptotic_func(x, d, e), 
                         start = list(e = 1),
                         control = nls.lm.control(maxiter = 1024))

bfit_asymptotic
asymptotic_coefs <- coef(bfit_asymptotic)

curve(asymptotic_func(x, d, asymptotic_coefs['e']), 
      add = TRUE, col = "gold1", lwd = 2)


##Gompertz fit
gompertz_func <- function(x, b, d, e) {
  1 + (d-1)* exp(-exp(-b*(x-e)))
}


bfit_gompertz <- nlsLM(y ~ gompertz_func(x, b, d, e), 
                       start = list(b = 0.1, e = 1),
                       control = nls.lm.control(maxiter = 1024))

bfit_gompertz
gompertz_coefs <- coef(bfit_gompertz)
gompertz_rss <- sum(resid(bfit_gompertz)^2)
gompertz_rss
curve(gompertz_func(x, gompertz_coefs['b'], d, gompertz_coefs['e']), 
      add = TRUE, col = "orange", lwd = 2)

##Modified Gompertz fit
mod_gompertz_func <- function(x, b, d, e) {
  1 + (d-1) * (1 - exp(-exp(b*(x-e))))
}
bfit_mod_gompertz <- nlsLM(y ~ mod_gompertz_func(x, b, d, e), 
                           start = list(b = 0.1, e = 1),
                           control = nls.lm.control(maxiter = 1024))

bfit_mod_gompertz
mod_gompertz_coefs <- coef(bfit_mod_gompertz)
mod_gompertz_rss <- sum(resid(bfit_mod_gompertz)^2)
mod_gompertz_rss
curve(mod_gompertz_func(x, mod_gompertz_coefs['b'], d, mod_gompertz_coefs['e']), 
      add = TRUE, col = "maroon", lwd = 2)


##Inverse sigmoidal fit
Inverse_sigmoidal_func <- function(x, d, b, e){((1-d) / (1 + exp(b*x-e)))+d}
bfit_inverse_sigm <- nlsLM(y ~ Inverse_sigmoidal_func(x, d, b, e), 
                           start = list(b = 0.1, e = 0.1),
                           control = nls.lm.control(maxiter = 1024))
bfit_inverse_sigm
inverse_sigm_coefs <- coef(bfit_inverse_sigm)
inverse_sigm_rss <- sum(resid(bfit_inverse_sigm)^2)
inverse_sigm_rss
curve(Inverse_sigmoidal_func(x, d, inverse_sigm_coefs['b'], inverse_sigm_coefs['e']), 
      add = TRUE, col = "blue", lwd = 2)


##Weibull I
weibull_I_func <- function(x, e, d, b) {
  1 + (d-1)*exp(-exp(-e*(log10(x)-log10(b))))
}

bfit_weibull_I <- nlsLM(y ~ weibull_I_func(x, e, d, b), 
                        start = list(e = 0.1, b = 0.1),
                        control = nls.lm.control(maxiter = 1024))
bfit_weibull_I
weibull_I_coefs <- coef(bfit_weibull_I)
weibull_I_rss <- sum(resid(bfit_weibull_I)^2)
weibull_I_rss
curve(weibull_I_func(x, weibull_I_coefs['e'], d, weibull_I_coefs['b']), 
      add = TRUE, col = "darkgreen", lwd = 2)


#Weibull II
weibull_II_func <- function(x, e, d, b) {
  1 + (d-1)*(1-exp(-exp(e*(log10(x)-log10(b)))))
}

bfit_weibull_II <- nlsLM(y ~ weibull_II_func(x, e, d, b), 
                         start = list(e = 2, b = 4.489996e+00),
                         control = nls.lm.control(maxiter = 1024))

bfit_weibull_II
weibull_II_coefs <- coef(bfit_weibull_II)
weibull_II_rss <- sum(resid(bfit_weibull_II)^2)

curve(weibull_II_func(x, weibull_II_coefs['e'], d, weibull_II_coefs['b']), 
      add = TRUE, col = "lightseagreen", lwd = 2)


##Lorentz I
d1 = max(y)
Lorentz_func <- function(x, e, d1, b) {
  1 + (d1-1)/(1+e*(x-b)^2) 
}

bfit_Lorentz <- nlsLM(y ~ Lorentz_func(x, e, d1, b), 
                      start = list(b = 1, e = 1),
                      control = nls.lm.control(maxiter = 1024))

bfit_Lorentz
Lorentz_coefs <- coef(bfit_Lorentz)
coef(bfit_Lorentz)
Lorentz_rss <- sum(resid(bfit_Lorentz)^2)
Lorentz_rss

curve(Lorentz_func(x, Lorentz_coefs['e'], d1, Lorentz_coefs['b']), 
      add = TRUE, col = "green2", lwd = 2)


##Lorentz II
d1 = max(y)
Lorentz_func_II <- function(x, e, d1, b) {
  d + (d1-d)/(1+e*(x-b)^2)
}


bfit_Lorentz_II <- nlsLM(y ~ Lorentz_func_II(x, e, d1, b), 
                      start = list(b = 1, e = 10),
                      control = nls.lm.control(maxiter = 1024))

bfit_Lorentz_II
Lorentz_coefs_II <- coef(bfit_Lorentz_II)
coef(bfit_Lorentz_II)
Lorentz_rss_II <- sum(resid(bfit_Lorentz_II)^2)
Lorentz_rss_II

curve(Lorentz_func_II(x, Lorentz_coefs_II['e'], d1, Lorentz_coefs_II['b']), 
      add = TRUE, col = "green4", lwd = 2)


##Bragg fit
d1 = max(y)
bragg_func <- function(x, e, b, d1) {
  1 + (d1-1)*exp(-e*(x-b)^2)
}

bfit_bragg <- nlsLM(y ~ bragg_func(x, e, d1, b), 
                    start = list(b = -1.599262e+04, e = 3.550939e-09),
                    control = nls.lm.control(maxiter = 1024))
bfit_bragg
bragg_coefs <- coef(bfit_bragg)
bragg_rss <- sum(resid(bfit_bragg)^2)
bragg_rss
curve(bragg_func(x, bragg_coefs['e'], d1, bragg_coefs['b']), 
      add = TRUE, col = "purple4", lwd = 2)



##Modified exponential decay
mod_exp_func <- function(x, e, b) {
  (d1 * exp(-e * x) + b * x)
}

bfit_mod_exp <- nlsLM(y ~ mod_exp_func(x, e, b), 
                        start = list(e = 0.1, b = 0.1), 
                        control = nls.lm.control(maxiter = 1024))

bfit_mod_exp
mod_exp_coefs <- coef(bfit_mod_exp)
mod_exp_rss <- sum(resid(bfit_mod_exp)^2)
mod_exp_rss
curve(mod_exp_func(x, mod_exp_coefs['e'], mod_exp_coefs['b']), 
      add = TRUE, col = "coral", lwd = 2)


##Gaussian-like peak with offset
Gaussian_Peak_Offset <- function(x, e, b) {
  d + (d1 - d) * exp(-e * (x - b)^2)
}

bfit_Gaussian_Peak <- nlsLM(y ~ Gaussian_Peak_Offset(x, e, b),
                            start = list(b = 0.1, e = 0.1),
                            control = nls.lm.control(maxiter = 1024))

##Extract coefficients
Gaussian_Peak_coefs <- coef(bfit_Gaussian_Peak)

##Calculate residual sum of squares
Gaussian_Peak_rss <- sum(resid(bfit_Gaussian_Peak)^2)
Gaussian_Peak_rss
curve(Gaussian_Peak_Offset(x, Gaussian_Peak_coefs['e'], Gaussian_Peak_coefs['b']),
      add = TRUE, col = "violet", lwd = 2)









##Plot the rss of linear fits versus rss of best-fitting models
str(combined_fits_all_df)


linear_func_data <- subset(combined_fits_all_df, Model_type == "Linear_func")

##Create a histogram with density plot
ggplot(linear_func_data, aes(x = Sum_of_squares + 0.000001)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "blue", alpha = 0.5) +
  geom_density(color = "red", size = 0.5) +
  theme_minimal() +
  scale_x_log10() +
  ylim(0, 3) +
  labs(title = "Sum of Squares for Linear Models",
       x = "Sum of Squares",
       y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))


##only subset for single lines per model, since many best-fitting models have 2 params 
ggplot(subset(best_fitting_models_df, Fit_param_name == "e"), aes(x = Sum_of_squares + 0.000001)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "blue", alpha = 0.5) +
  geom_density(color = "red", size = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  #ylim(0, 45) +
  ylim(0, 3) +
  labs(title = "Sum of Squares for Best-fitting Models",
       x = "Sum of Squares",
       y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))



##Non-logged:
ggplot(linear_func_data, aes(x = Sum_of_squares)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "blue", alpha = 0.5) +
  #geom_density(color = "red", size = 0.5) +
  theme_minimal() +
  ylim(0, 45) +
  labs(title = "Sum of Squares for Linear Models",
       x = "Sum of Squares",
       y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))


##only subset for single lines per model, since many best-fitting models have 2 params 
ggplot(subset(best_fitting_models_df, Fit_param_name == "e"), aes(x = Sum_of_squares)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "blue", alpha = 0.5) +
  #geom_density(color = "red", size = 0.5) +
  theme_minimal() +
  ylim(0, 45) +
  labs(title = "Sum of Squares for Best-fitting Models",
       x = "Sum of Squares",
       y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))



#As used in the Thesis:
##Overlapping histograms emphasising benefit of using non-linear models to account for drug effect as a function of concentration:
# Add an identifier to each dataset
linear_func_data$type <- "Linear models"
best_fitting_data <- subset(best_fitting_models_df, Fit_param_name == "e")
best_fitting_data$type <- "Best-fitting Models"

##Combine the two datasets
combined_data <- rbind(
  linear_func_data[, c("Sum_of_squares", "type")],
  best_fitting_data[, c("Sum_of_squares", "type")]
)

##overlapping histogram
ggplot(combined_data, aes(x = Sum_of_squares + 0.000001, fill = type)) +
  geom_histogram(aes(y = ..density..), bins = 100, alpha = 0.5, position = "identity") +
  geom_density(aes(color = type), size = 0.5, alpha = 0.3) +
  scale_x_log10() +
  ylim(0, 0.8) +
  theme_minimal() +
  labs(title = "",
       x = "Log10 (Sum of Squares +1e-06)",
       y = "Density",
       fill = "Model Type",
       color = "Model Type") +
  theme(plot.title = element_text(hjust = 0.5))




#As used in Thesis:
### Best-fitting models by Species and by drugs
best_fitting_models_df2 <- subset(best_fitting_models_df1, Fit_param_name == "e" & Species != "P_copri" & Species != "B_longum" & Species != "B_wadsworthia")
species_model_table <- table(best_fitting_models_df2$Species, best_fitting_models_df2$Model_type)
drug_model_table <- table(best_fitting_models_df2$drug, best_fitting_models_df2$Model_type)


#plotting the count of fits across model types:
drug_df <- as.data.frame(drug_model_table)

##Rename the columns for clarity (Var2 represents Model_type)
colnames(drug_df) <- c("Drug", "Model_type", "Count")

##Aggregate counts by Model_type
drug_summary <- aggregate(Count ~ Model_type, data = drug_df, sum)
sum(drug_summary$Count) #this may not sum to 30*21 entirely because some models return NA with fitting with set/defaulted starting vals for 'e' and 'b' params 

drug_summary$Model_type <- reorder(drug_summary$Model_type, -1*(drug_summary$Count))

ggplot(drug_summary, aes(x = Model_type, y = Count)) +
  geom_bar(stat = "identity", fill = "grey50") +
  labs(title = "Counts of best model fits by model type",
       x = "Model type",
       y = "Total Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




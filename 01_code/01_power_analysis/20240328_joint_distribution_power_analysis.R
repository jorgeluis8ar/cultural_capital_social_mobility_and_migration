# Description ------------------------------------------------------------------

#'
#'Script: 20240312_conjoint_analysis_joint_distribution_power_calc
#'
#' RA1: Jorge Luis Ochoa Rincon
#' RA2: Hernan Carvajal
#' 
#' Date created: 2024-03-12
#' Date modified: 2024-03-12
#' Objective:
#' Demo and power calculation for the conjoint analysis using both a uniform
#' and a joint distribution. We follow de la Cuesta, B., et al (2021) paper and the
#' code provided in the package factorEx to produce our results. The following are
#' of the used material:
#' Paper: de la Cuesta, B., et al (2021)- Improving the External Validity of Conjoint Analysis: The Essential Role of Profile Distribution
#' Code demo: https://github.com/naoki-egami/factorEx
#' EGAP demo: https://egap.org/resource/script-power-analysis-simulations-in-r/
#' 
#' Notes
#' We do not use special characters in the code, except for values in the data set
#' in order to harmonize output in all operative systems

# ---------------------------------------------------------------------------- #

# 0. Cleaning the console and environment --------------------------------------

cat("\f")                           # Cleaning the console
rm(list=ls())                       # Cleaning the environment
options("scipen"=100, "digits"=4)   # Format and number of digits for number
Sys.setlocale("LC_TIME", "C")       # English format for dates

# 1. Packages ------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
devtools::install_github("m-freitag/cjpowR")                            # Installing the package from source
devtools::install_github("naoki-egami/factorEx", dependencies=TRUE)     # Installing the package from source
pacman::p_load(tidyverse,tidytext,topicmodels,SnowballC,dplyr,tidyr,fixest, tictoc,
               devtools,cjpowR,ggplot2,factorEx,randomizr,data.table,rcartocolor)

clrs <- carto_pal(name = "Prism")

# 2. Directory -----------------------------------------------------------------

if (Sys.info()["user"] == "jl.ochoa") {
  
  conjoint_folder <- "C:/Users/jl.ochoa/OneDrive - Universidad de los Andes/Proyecto Hospitales/"
  
}else if (Sys.info()["user"] == "username") {
  
  conjoint_folder <- "Change directory to the Dropbox/Onedrive/Drive folder"
  
}

# joint distributions 
joint <- fread(paste0(conjoint_folder,"01_datos/raw/joint_dis.csv"),encoding = "UTF-8")

# Simulation -------------------------------------------------------------------

#### Creating the partial probabilites - joint and marginal ####

# Marginal uniform distribution for the controls
experience <- c("Baja" = 1/3, 
                "Promedio" = 1/3, 
                "Alta" = 1/3)
big.five <- c("Extrovertido(a) y Sociable" = 1/5, 
              "Amable y Generoso(a)" = 1/5,
              "Minucioso(a) y Eficiente" = 1/5,
              "Calmado(a)" = 1/5,
              "Ingenioso(a) y Analítico(a)" = 1/5)
region <- c("Bogotá" = 1/2,
            "Santiago" = 1/2)
clase <- c("Baja" = 1/2,
           "Alta" = 1/2)

#  Joint distribution for gender, income and education
table.joint <- joint %>% 
  summarise(prob = sum(perc,na.rm = T),.by = c("sex_num","inc_num","educ_num")) %>%
  mutate(prob = (prob/100)/2) %>% 
  rename("income" = "inc_num",
         "education" = "educ_num",
         "gender" = "sex_num") %>% 
  xtabs(prob ~ income + education + gender, .)
table.joint <- table.joint[,unique(joint$educ_num),]
table.joint <- table.joint[unique(joint$inc_num),,]
table.joint <- table.joint[,,unique(joint$sex_num)]

# Final object with the target distribution for the simulation
target_dist_simulation2 <- list(`income:education:gender` = ftable(table.joint) %>% as.table(),
                                char_experience = experience,
                                char_big_five = big.five,
                                char_region = region,
                                char_clase = clase)
target_dist_simulation2


partial_joint_name_test2  <- list(c("income", "education", "gender"), 
                                  "char_experience","char_big_five", "char_region","char_clase")

#### Creating the structure for the simulation ####

#' Pool of randomly drawn attributes. All included here are the controls with 
#' marginal uniform distributions. This includes:
#' 
#' 1 Experience:
#' 2. Region
#' 3. Clase:
#' 4. Big five
#' 
df.random <- expand.grid(c("Baja","Promedio","Alta"),
                         c("Bogotá","Santiago"),
                         c("Baja","Alta"),
                         c("Extrovertido(a) y Sociable", "Amable y Generoso(a)","Minucioso(a) y Eficiente","Calmado(a)","Ingenioso(a) y Analítico(a)"))
df.random <- as.data.frame(df.random)
colnames(df.random) <- c("char_experience","char_region","char_clase","char_big_five")

#' Now, we have to sample attributes accordingly to the joint distribution provided.
#' This, I use the proportion of each combination to match to each level of N size.

joint$num_prop <- (joint$perc/2)/100
# Get the number of repetitions for each combination

list.amce.power.sample <- list()
sample.indicator <- 1



tic()
for (n.conjoint in seq(1500,2000,1)) {
  
  total_size <- n.conjoint*2*6
  
  cat("\f")
  cat("\nSample size:", total_size)

  df.sim <- joint %>% 
    mutate(samples = round(total_size*num_prop)) %>% 
    select(sex_num, educ_num, inc_num,samples)
  
  # Select and repeat the number of attributes combinations
  df.sim <- as.data.frame(lapply(df.sim, rep, df.sim$samples))
  df.sim <- df.sim %>% select(sex_num, educ_num, inc_num)
  
  # if number are not exactly the same, randomly add/remove observations
  if (nrow(df.sim) < total_size) {
    
    cat("\nLess observations\n")
    
    obs.missing <- total_size-nrow(df.sim)
    df.remaining <- joint %>% 
      mutate(ran_number = runif(n = nrow(joint),min = 0,max = 1)) %>% 
      select(sex_num, educ_num, inc_num,ran_number) %>% 
      arrange(desc(ran_number)) %>% 
      slice(1:obs.missing) %>% 
      select(sex_num, educ_num, inc_num)
    
    df.sim <- dplyr::bind_rows(df.sim,df.remaining)
  }else if (nrow(df.sim) > total_size) {
    cat("\nMore observations\n")
    df.sim <- df.sim[sample(nrow(df.sim), total_size),]
  }  
  
  # Pasting the marginal uniform attributives
  df.sim$random_obs <- round(runif(n = nrow(df.sim),min = 1,max = nrow(df.random)))
  df.sim <- dplyr::bind_cols(df.sim,df.random[df.sim$random_obs,])
  df.sim <- df.sim %>% select(sex_num, educ_num, inc_num,char_experience,char_region,char_clase,char_big_five)
  df.sim <- df.sim %>% 
    rename("gender" = "sex_num",
           "education" = "educ_num",
           "income" = "inc_num")
  df.sim$gender <- factor(df.sim$gender,levels = unique(joint$sex_num))
  df.sim$education <- factor(df.sim$education,levels = unique(joint$educ_num))
  df.sim$income <- factor(df.sim$income,levels = unique(joint$inc_num))
  
  df.sim <- as.data.table(df.sim)
  
  number_reps = nrow(df.sim)/12
  
  df.sim[,id := rep(seq(1,n.conjoint,1),12)]
  setorder(df.sim,id)
  df.sim$pair_id <- rep(seq(1,nrow(df.sim)/2,1),2) %>% sort(decreasing = FALSE)
  setorder(df.sim,id,pair_id)
  
  df.sim <- as.data.frame(df.sim)
  
  list.amce.power <- list()
  amce.indicator <- 1
  for (amce in seq(0.1, 0.1, length.out = 1)) {
    
    list.reg <- parallel::mclapply(seq(1,10000,1), function(x){
      
      cat("\nAMCE:", amce)  
      cat("\nSimulation:", x, "with sample size:",n.conjoint,"\n")
      
      df.sim$outcome <- rbinom(n=nrow(df.sim), size=1, prob=.5)          # Do a random assignment
      
      # Now, before running the simulation, we gotta make sure that we will get the desired amce
      n_replace <- round(amce*n.conjoint*6)
      cat("\nObservations to replace:",n_replace,"\n")
      
      n_alta <- nrow(df.sim[(df.sim$outcome == 1)&(df.sim$char_clase == "Alta"),])
      n_baja <- nrow(df.sim[(df.sim$outcome == 1)&(df.sim$char_clase == "Baja"),])
      
      n_desired <- n_baja + n_replace
      n_2_replace <- n_desired - n_alta
      n_2_replace <- n_2_replace*runif(n = 1,min = 0,max = 1)
      n_2_replace <- round(n_2_replace)
      
      if (n_2_replace < 0) {
        
        index_2_replace <- sample(grep("1",x = df.sim$outcome)[grep("1",x = df.sim$outcome) %in% grep("Alta",x = df.sim$char_clase)],abs(n_2_replace))
        df.sim$outcome2 <- df.sim$outcome
        df.sim$outcome2[index_2_replace] <- 0  
      }else if (n_2_replace > 0) {
        
        index_2_replace <- sample(grep("0",x = df.sim$outcome)[grep("0",x = df.sim$outcome) %in% grep("Alta",x = df.sim$char_clase)],n_2_replace)
        df.sim$outcome2 <- df.sim$outcome
        df.sim$outcome2[index_2_replace] <- 1
      }else if (n_2_replace == 0) {
        df.sim$outcome2 <- df.sim$outcome  
      } 
      
      # Running the command    
      sim_pamce <- design_pAMCE(formula            = outcome2 ~ char_clase + income + gender + education  + char_experience + char_region + char_big_five,
                                factor_name        = c("char_clase"),
                                data               = df.sim,
                                pair_id            = df.sim$pair_id,
                                cluster_id         = df.sim$id,
                                target_dist        = target_dist_simulation2, 
                                target_type        = "partial_joint",
                                partial_joint_name = partial_joint_name_test2)
      
      # Saving the data frame with the results
      df.pamce <- sim_pamce$AMCE$char_clase
      
      # Creating the z-value
      acento_pamce <- df.pamce %>% filter(type == "target") %>% pull(estimate)
      acento_pamce_se <- df.pamce %>% filter(type == "target") %>% pull(se)
      
      # Getting the associated p value
      z_score_amce <- acento_pamce/acento_pamce_se
      p_value_amce <- 2*pnorm(q=z_score_amce, lower.tail=FALSE)
      power <- pnorm(z_score_amce - qnorm(1 - alpha/2)) + pnorm(-z_score_amce - qnorm(1 - alpha/2))
      
      # Saving results in the data.frame
      df.pamce <- as.data.table(df.pamce)
      df.pamce <- df.pamce[type == "target"]
      df.pamce[, z_score := z_score_amce]
      df.pamce[, p_value := p_value_amce]
      df.pamce[, ind_pvalue := ifelse(p_value <= 0.05,TRUE,FALSE)]
      df.pamce[, power_pamce := power]
      return(df.pamce)
    })
    list.reg <- plyr::rbind.fill(list.reg %>% data.table::rbindlist(use.names = T,fill = T))
    
    amce.df <- data.frame(amce = amce,
                          power = mean(list.reg$ind_pvalue,na.rm = T),
                          mean.coef = mean(list.reg$estimate,na.rm = T),
                          mean.power = mean(list.reg$power_pamce,na.rm = T))
    
    list.amce.power[[amce.indicator]] <- amce.df
    amce.indicator <- amce.indicator + 1
  }
  
  list.amce.power <- plyr::rbind.fill(list.amce.power %>% data.table::rbindlist(use.names = T,fill = T))
  list.amce.power$sample_size <- n.conjoint
  
  list.amce.power.sample[[sample.indicator]] <- list.amce.power
  sample.indicator <- sample.indicator + 1
  saveRDS(object = list.amce.power.sample,file = paste0("c:/Users/jl.ochoa/OneDrive - Universidad de los Andes/pamce2/pamce_simulation_samples_fith_",n.conjoint,".rds"))
}
time <- toc()
list.amce.power.sample <- plyr::rbind.fill(list.amce.power.sample %>% data.table::rbindlist(use.names = T,fill = T))
# saveRDS(object = list.amce.power.sample,file = "c:/Users/jl.ochoa/OneDrive - Universidad de los Andes/pamce_simulation_2000_first_20.rds")

ggplot(data = list.amce.power.sample) +
  geom_line(aes(x = amce, y = power, group = factor(sample_size), color = factor(sample_size), linetype = factor(sample_size))) +
  scale_x_continuous(breaks = seq(0, 0.2, 0.025)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  geom_hline(yintercept=0.8, linetype="dashed", color = clrs[8]) +
  scale_linetype_manual(values = c("1500" = "solid","1625" = "solid",
                                   "1750" = "solid", "1875" =  "solid",
                                   "2000" = "solid", "solid")) +
  scale_color_manual(values = c("1500" = clrs[2],"1625" = clrs[5],
                                "1750" = clrs[9],"1875" = clrs[6],
                                "2000" =  clrs[12], "blue")) + 
  labs(x = "AMCE", y = "Power", color = "Sample size", linetype = "Sample size") + 
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 10),axis.title = element_text(color = "black", size = 14),
        legend.position = c(0.1,0.68), legend.background = element_blank(), legend.key = element_blank(), 
        legend.text = element_text(color = "black", size = 12), legend.title = element_text(color = "black", size = 15),
        panel.grid.minor = element_blank())

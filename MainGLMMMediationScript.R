#################################################
# CODE TO CONDUCT MEDIATION AND OTHER ANALYSES  #
# EFFECT OF RODENTS ON COCCIDIOIDES AS MEDIATED #
# BY BURROWS, CARRIZO PLAIN, CA                 #
# CREATED BY: JENNIFER HEAD                     #
# LAST EDITED: FEB 15, 2023                     #
#################################################

####################################
## PART 0. LOAD DATA AND PACKAGES ##
####################################

package_list <-c("ggplot2", "dplyr", "paletteer", "mice", "lme4", 
                 "stargazer", "table1", "stringr", "doParallel", "foreach", 
                 "here")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(package_list, character.only=TRUE)

# Load the data
setwd(here())
dat <- read.csv("cleaned_data_plus_results_Mar27_2023.csv")

#################################
## PART I. DESCRIPTIVE RESULTS ##
#################################

# Overall prevalence
(tot <- sum(!is.na(dat$positive))) # total samples analyzed
(pos <- sum(dat$positive == 1, na.rm = T)) # total samples positive
(posrate <- pos/tot) # percent positive

# Make some data tables
prev_by_trtmnt <- dat %>% group_by(factorial) %>% summarize(n = sum(!is.na(positive)),
                                                     count = sum(positive == 1, na.rm = T), 
                                                     percent = mean(positive, na.rm = T))

prev_by_season <- dat %>% group_by(season) %>% summarize(n = sum(!is.na(positive)),
                                                   count = sum(positive == 1, na.rm = T), 
                                                   percent = mean(positive, na.rm = T))

prev_by_trt_ssn <- dat %>% group_by(factorial, season) %>% summarize(n = sum(!is.na(positive)),
                                                              count = sum(positive == 1, na.rm = T), 
                                                              percent = mean(positive, na.rm = T))


## Make a table 1 style table using the table1 package
table1 <- table1( ~ factorial | season,
                  data = dat)

table1

## Make a barplot of overall positive rates
ggplot(prev_by_trtmnt) + 
  geom_bar(aes(x = factorial, y = percent*100), stat = "identity", fill = "darkcyan") +
  theme_bw() + xlab("") + ylab("Percent of samples with Coccidioides") +
  theme(text = element_text(size = 18))

#ggsave("PercentSamples.jpg", width = 6, height = 5)

## Make a barplot of overall positive rates by season
ggplot(prev_by_trt_ssn) + 
  geom_bar(aes(x = factorial, y = percent*100, fill = season), 
                         stat = "identity", position = "dodge") +
  theme_bw() + xlab("") + ylab("Percent of samples with Coccidioides") +
  theme(text = element_text(size = 18)) +
  scale_fill_paletteer_d(`"awtools::a_palette"`)


## Supplemental figure
prev_by_trt_ssn$season2 <- factor(prev_by_trt_ssn$season,
                        levels = c("Spring", "Summer", "Fall", "April"),
                        labels = c("Spring '21", "Summer '21", "Fall '21", "Spring '22"))

prev_by_trt_ssn$factorial <- factor(prev_by_trt_ssn$factorial,
                          levels = c("SE", "SN", "PE", "PN"),
                          labels = c("Surface without rodents",
                                     "Surface with rodents",
                                     "Precincts without rodents",
                                     "Precincts with rodents"))

## Make a barplot of overall positive rates by season
ggplot(prev_by_trt_ssn) + 
  geom_bar(aes(x = factorial, y = percent*100, fill = season2), 
                         stat = "identity", position = "dodge") +
  theme_bw() + xlab("") + ylab("Percent of samples with Coccidioides") +
  theme(text = element_text(size = 18)) +
  scale_fill_paletteer_d(`"awtools::a_palette"`, name = "Season") +
  facet_wrap(~season2) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))

#ggsave("PercentSamples_bySeason.jpg", width = 12, height = 7)

#######################################################
## PART II. MULTIPLE IMPUTATION OF MISSING VARIABLES ##
#######################################################

### Examine variables
summary(dat$positive) # missing outcome
summary(dat$vwc_med) # missing water content (if soil too hard for probe)

#### Do multiple imputation for the vwc ####
missing <- dat %>% dplyr::select(c("vwc", "plot_num", "treatment", "exclosure", 
                                      "veg_level", "soil_temp", "positive"))
summary(missing)

#Do the imputation
vwc_imp <- mice(missing,
                m = 5)

#Examine what has been imputed
vwc_imp$imp$vwc
vwc_imp$imp$positive

#Examine 1 of 5 full sets
missing1 <- complete(vwc_imp,1)
missing2 <- complete(vwc_imp,2)
missing3 <- complete(vwc_imp,3)
missing4 <- complete(vwc_imp,4)
missing5 <- complete(vwc_imp,5)

#Compute mean of imputed values
missing_mat <- cbind(missing1$vwc, missing2$vwc, missing3$vwc, missing4$vwc, missing5$vwc)
mean_vwc_imp <- matrixStats::rowMedians(missing_mat)
dat$vwc_imp <- mean_vwc_imp 

mean(dat$vwc_imp)
sd(dat$vwc_imp)
range(dat$vwc_imp)

mean(dat$soil_temp)
sd(dat$soil_temp)
range(dat$soil_temp)
tail(sort(dat$soil_temp))

mean(dat$veg_level, na.rm = T)
sd(dat$veg_level, na.rm = T)
range(dat$veg_level, na.rm = T)

#Compare to observed
dat %>% group_by(treatment, exclosure) %>% summarize(mean(vwc,na.rm = T))
dat %>% group_by(treatment, exclosure) %>% summarize(mean(vwc_imp,na.rm = T))
missing1 %>% group_by(treatment, exclosure) %>% summarize(mean(vwc,na.rm = T))

# Supplemental figure: boxplot of soil moisture

dat$season2 <- factor(dat$season,
                      levels = c("Spring", "Summer", "Fall", "April"),
                      labels = c("Spring '21", "Summer '21", "Fall '21", "Spring '22"))


exc.labs <- c("Exclosure", "Non-exclosure")
names(exc.labs) <- c("E", "N")

ggplot(dat) + 
  geom_boxplot(aes(x = season2, y = vwc_imp, fill = treatmentX)) + 
  facet_wrap(~exclosureX, labeller = labeller(exclosureX = exc.labs)) +
  xlab("") + ylab("Volumetric water content (%)") + theme_bw() +
  theme(legend.position = "top") +
  scale_fill_brewer("Burrow", palette = "Paired", labels = c("Yes", "No"))

#ggsave(file = "MoistureComparison_all_sns.jpg", dpi = 250, width = 8, height = 5)

#####################################################
## PART IIIA. CREATE GLMMS FOR RODENTS AND BURROWS ##
#####################################################

# relevel
dat$factorial <- factor(dat$factorial, levels = c("SE", "SN", "PE", "PN"))

# Set April to season (dat)
dat$season2 <- dat$season
dat$season2[dat$season == "April"] <- "Spring"
dat$season <- factor(dat$season, levels = c("Spring", "Summer", "Fall", "April"))

##### Establish variables for sensitivity analyses -- edit down below #####

# Sensitivity analyses 1 - which exclosures may be breached #
dat$active_exclosure_repeats <- as.numeric(dat$plot %in% c("C10", "S07"))
dat$active_exclosure_any <- as.numeric(dat$plot %in% c("C10", "S07", "C08", "C04", "S04"))

# Sensitivity analysis 2 - recode positive #
dat$positive_1of4 <- as.numeric(dat$Ct_detects >= 1)
dat$positive_2of4 <- as.numeric(dat$Ct_detects >= 2)
dat$positive_4of4 <- as.numeric(dat$Ct_detects >= 4)

##############################################
#### EDIT HERE FOR SENSITIVITY ANALYSES ######
dat_orig <- dat

## Options for sensitivity analysis - comment out for main analysis
# dat <- dat %>% subset(active_exclosure_repeats != 1)
# dat <- dat %>% subset(active_exclosure_any != 1)
# dat$positive <- dat$positive_1of4
# dat$positive <- dat$positive_2of4
# dat$positive <- dat$positive_4of4

dat <- dat_orig ## reset

##############################################

##############################################

# Make a GLMM, adjusting for various factors
mod_unadj_noint <-  glmer(positive ~ season + 
                            exclosure + treatment + (1|plot_num/cluster_pair), 
                          data = dat, family = "binomial")

summary(mod_unadj_noint)

mod_adj_noint <- glmer(positive ~ vwc_imp + season +
                         exclosure + treatment + (1|plot_num/cluster_pair), 
                       data = dat, family = "binomial")

summary(mod_adj_noint)

# Make a nice table of results -- it will save as html output
stargazer(mod_unadj_noint, mod_adj_noint, align = T, ci = T, apply.coef = exp, p.auto = F,
          order = c(5,6,2,3,4,1,7),
          column.labels = c("Model A", "Model B"),
          covariate.labels = c("Rodents present (ref: absent)", 
                               "Burrow (ref: surface)", 
                               "Summer 2021 (ref: dat 2021)",
                               "Fall 2021 (ref: dat 2021)",
                               "dat 2022 (ref: dat 2021)",
                               "Volumetric water content (%)", 
                               "Intercept"),
          dep.var.labels = "OR (95% CI)",
          dep.var.caption = "Presence of \\textit{Coccidioides} in soils",
          out="models.htm")

#####################################################
## PART IIIB. CREATE GLMMS FOR SOIL CONDITIONS     ##
#####################################################

# soil moisture
mod_vwc <-lmer(vwc_imp ~ season + 
                          exclosure + treatment + (1|plot_num/cluster_pair), 
                          data = dat)

summary(mod_vwc)
confint(mod_vwc)

# ordinal soil level
mod_veg <-lmer(veg_level ~ season + 
                   exclosure + treatment + (1|plot_num/cluster_pair), 
                 data = dat)

summary(mod_veg)
confint(mod_veg)

# soil temperature
mod_tmp <- lmer(soil_temp ~ season + 
                   exclosure + treatment + (1|plot_num/cluster_pair), 
                 data = dat)

summary(mod_tmp)
confint(mod_tmp)

#################################################
##  PART IVA. MEDIATION ANALYSIS               ##
##      EXPOSURE (X.A) : RODENTS               ##
##      MEDIATOR (M) : BURROWS                 ##
##      OUTCOME  (Y) : COCCIDIOIDES POS        ##
##  Approach following G-computation approach  ##
##  See Arah and Wang, 2015 for more details   ##
#################################################

# Step 0. Define key variables
# A is the main effect (rodents, exclosure == 1)
# M is the mediator (burrows, treatment == 1)
# Y is the outcome (Coccidioides, positive == 1)

dat$A <- dat$exclosure
dat$M <- dat$treatment
dat$Y <- dat$positive
dat_full <- dat %>% subset(!is.na(Y))

#Any subsetting by season done here -- comment out for full analysis
# dat_full <- dat_full %>% subset(season == "April")

#Step 1a: obtain appropriate distribution of each variable
P.X <- mean(dat_full$A)

#Step 1b: Model the mediator on exposure and relevant covariates
alpha0 <- 0 #when A is 0, M is 0
alpha1 <- 0.5 # when A is 1, M is 0.6 (60% chance of being a 1)
RMSE.M <- 0.1

#Step 1c: Chose which model to run! Don't run them both!

# Model 1. Model outcome on A, M 
E.Y <- glmer(Y ~ season + A + M + (1|plot_num/cluster_pair), data = dat_full, family = "binomial")

# Model 2 -- including other variables
E.Y <- glmer(Y ~ vwc_imp + season + A + M + (1|plot_num/cluster_pair), 
             data = dat_full, family = "binomial")

coefs <- fixef(E.Y)

#Step 2a. Create J copies of the original sample

ncopies <- 1
Copy <- dat_full

for (rep in 2:ncopies){
  Copy <- rbind(Copy, dat_full) 
}

plots <- unique(paste0(dat_full$plot)) # get unique Carrizo plots
nboot <- 1000 # number of bootstrap iterations

# a function to predict Y's
get_y <- function(Xvar,Mvar, binom = T){ 
  
  # Recreate dataset for prediction
  pred_df <- boot.df[,c("vwc_imp", "season", Xvar, Mvar)]
  colnames(pred_df) <- c("vwc_imp", "season", "A", "M")
  
  # First, we predict the mean of Y, then we get the SD, then we resample, accounting for the SD
  # Predict the mean of Y
  pred_df$Y <- predict(E.Y, pred_df, re.form = NA)
  
  # Get the SD 
  mm <- model.matrix(terms(E.Y),pred_df)
  pvar1 <- sqrt(diag(mm %*% tcrossprod(vcov(E.Y),mm)))
  
  # Sample from normal distribution, with mean as given and sd as found
  vals <- pmin(rnorm(n, mean = pred_df$Y, sd = pvar1), 0) ##??why are some > 0?
  
  # Sample from binomial distribution with prob of success equal to the mean probability
  if (binom == T){
    vals <- rbinom(n = n, size = 1, prob = exp(vals))
  }
  
  return(vals)
}

# set up parallel computing
numCore <- detectCores() - 1
cl <- makeCluster(numCore)
registerDoParallel(cl)

# parallel loop for the bootstrapping
res <- foreach (i = 1:nboot, .combine = rbind,
                .packages = c('dplyr', 'lme4')) %dopar% {
                  
                  # create a boostrapped dataframe, boot.df
                  boot.df <- dat_full %>% subset(positive == 20) #empty dataframe with same structure as spring_full
                  
                  # Bootstrapped samples must obey clustering - so drawing based on plots
                  
                  # Sample the plots 20 times, with replacement
                  boot_samp <- sample(plots, 20*ncopies, replace = T)
                  plots_samp <- sort(unique(boot_samp)) # the plots sampled
                  boot_table <- table(boot_samp) # a table of how often each plot was sampled
                  pn <- length(plots_samp) # the number of unique plots sampled
                  
                  # A loop to recreated the dataframe, by appending the dataset sampled plot at a time
                  for (p in 1:pn){
                    boot_df_i <- subset(Copy, plotBC == plots[p])
                    plot_n <- boot_table[p]
                    
                    k <- 1
                    while (plot_n >= k){
                      boot.df <- rbind(boot_df_i, boot.df)
                      k <- k+1
                    }
                  }
                  
                  n <- nrow(boot.df) # number of rows in bootstrapped dataframe
                  
                  #Step 2b. Simulate an intervention variable that is marginally independent of the simulated covariates
                  boot.df$X <- rbinom(n, 1, P.X)
                  boot.df$X1 <- 1
                  boot.df$X0 <- 0
                  
                  #Step 2c. Simulate each M as a function of its parents
                  P.Mx <- alpha0 + alpha1*boot.df$X + rnorm(n, mean = 0, sd = RMSE.M*boot.df$X)
                  P.M1 <- alpha0 + alpha1 + rnorm(n, mean = 0, sd = RMSE.M)
                  boot.df$Mx <- rbinom(n, 1, P.Mx)
                  boot.df$M1 <- rbinom(n, 1, P.M1)
                  boot.df$M0 <- rbinom(n, 1, alpha0)
                  
                  #Step 2d. simulate a Y for each type of effect as a function of its parents
                  
                  boot.df$Y_TE <- get_y("X", "Mx", binom = T)
                  boot.df$Y_PDE <- get_y("X", "M0", binom = T)
                  boot.df$Y_TIE <- get_y("X1", "Mx", binom = T)
                  boot.df$Y_TDE <- get_y("X", "M1", binom = T)
                  boot.df$Y_PIE <- get_y("X0", "Mx", binom = T)
                  
                  
                  #Step 3. Regress each potential health outcome against X
                  TE  <- coef(glm(Y_TE  ~ X, data = boot.df, family = "binomial"))[2]
                  PDE <- coef(glm(Y_PDE ~ X, data = boot.df, family = "binomial"))[2]
                  TIE <- coef(glm(Y_TIE ~ X, data = boot.df, family = "binomial"))[2]
                  TDE <- coef(glm(Y_TDE ~ X, data = boot.df, family = "binomial"))[2]
                  PIE <- coef(glm(Y_PIE ~ X, data = boot.df, family = "binomial"))[2]
                  
                  #Print summary
                  res <- c("TE" = TE, "PDE" = PDE, "TIE" = TIE, "TDE" = TDE, "PIE" = PIE)
                  
                  res
                }

stopCluster(cl)

## Summarize the results - median and 95% CI
MA_res <- matrixStats::colQuantiles(exp(res), probs = c(0.025, 0.5, 0.975), na.rm = T)
rownames(MA_res) <- c("TE", "PDE", "TIE", "TDE", "PIE")
MA_res

# percent mediated by burrows = TIE/TE
quantile(res[,3]/(res[,3] + res[,4])*100, probs = c(0.025, 0.5, 0.975), na.rm = T)
100-quantile(res[,3]/(res[,3] + res[,4])*100, probs = c(0.025, 0.5, 0.975), na.rm = T)


#################################################
##  PART IVB. MEDIATION ANALYSIS               ##
##      EXPOSURE (X/A) : BURROWS               ##
##      MEDIATOR (M) : SOIL CHARACTERISTICS    ##
##      OUTCOME  (Y) : COCCIDIOIDES POS        ##
##  Approach following G-computation approach  ##
##  See Arah and Wang, 2015 for more details   ##
#################################################

# Step 0. Define key variables
# A is the main effect (burrows, treatment == 1)
# M is the mediator (soil conditions, i.e., vwc_imp)
# Y is the outcome (Coccidioides, positive == 1)

dat$A <- dat$treatment
dat$M <- dat$vwc_imp
dat$Y <- dat$positive
dat_full <- dat %>% subset(!is.na(Y))

#Any subsetting by season done here -- comment out for full analysis
# dat_full <- dat_full %>% subset(season == "April")

#Step 1a: obtain appropriate distribution of each variable
P.X <- mean(dat_full$A)

#Step 1b: Model the mediator on exposure and relevant covariates
E.M <- lmer(vwc_imp ~ A + 
                   season + exclosure + 
                   (1|cluster_pair), 
                data = dat)

coefsM <- fixef(E.M)

#Step 1c: Model Y dependent on A and M

# Model 2 -- including other variables
E.Y <- glmer(Y ~ exclosure + season + A + M + (1|plot_num/cluster_pair), 
             data = dat_full, family = "binomial")

coefsY <- fixef(E.Y)

#Step 2a. Create J copies of the original sample

ncopies <- 1
Copy <- dat_full

for (rep in 2:ncopies){
  Copy <- rbind(Copy, dat_full) 
}

plots <- unique(paste0(dat_full$plot)) # get unique Carrizo plots
nboot <- 1000 # number of bootstrap iterations

# a function to predict Y's -- do need to repeat from above
get_y <- function(Xvar,Mvar, binom = T){ 

  # Recreate dataset for prediction
  pred_df <- boot.df[,c("exclosure", "season", Xvar, Mvar)]
  colnames(pred_df) <- c("exclosure", "season", "A", "M")
  
  # First, we predict the mean of Y, then we get the SD, then we resample, accounting for the SD
  # Predict the mean of Y
  pred_df$Y <- predict(E.Y, pred_df, re.form = NA)
  
  # Get the SD 
  mm <- model.matrix(terms(E.Y),pred_df)
  pvar1 <- sqrt(diag(mm %*% tcrossprod(vcov(E.Y),mm)))
  
  # Sample from normal distribution, with mean as given and sd as found
  vals <- pmin(rnorm(n, mean = pred_df$Y, sd = pvar1), 0) ##??why are some > 0?
  
  # Sample from binomial distribution with prob of success equal to the mean probability
  if (binom == T){
    vals <- rbinom(n = n, size = 1, prob = exp(vals))
  }
  
  return(vals)
}

# A function to get M's
get_m <- function(Xvar, binom = F){ 
  
  # Recreate dataset for prediction
  pred_df <- boot.df[,c("exclosure", "season", Xvar)]
  colnames(pred_df) <- c("exclosure", "season", "A")
  
  # First, we predict the mean of Y, then we get the SD, then we resample, accounting for the SD
  # Predict the mean of Y
  pred_df$M <- predict(E.M, pred_df, re.form = NA)
  
  # Get the SD 

  pvar1 <- sqrt(var(as.vector(fixef(E.M) %*% t(getME(E.M,"X")))))

  # Sample from normal distribution, with mean as given and sd as found
  vals <- pmax(rnorm(n, mean = pred_df$M, sd = pvar1), 0) 
  
  # Sample from binomial distribution with prob of success equal to the mean probability
  if (binom == T){
    vals <- rbinom(n = n, size = 1, prob = exp(vals))
  }
  
  return(vals)
}

# set up parallel computing
numCore <- detectCores() - 1
cl <- makeCluster(numCore)
registerDoParallel(cl)

# parallel loop for the bootstrapping
res <- foreach (i = 1:nboot, .combine = rbind,
                .packages = c('dplyr', 'lme4')) %dopar% {
                
                  # create a boostrapped dataframe, boot.df
                  boot.df <- dat_full %>% subset(positive == 20) #empty dataframe with same structure as spring_full
                  
                  # Bootstrapped samples must obey clustering - so drawing based on plots
                  
                  # Sample the plots 20 times, with replacement
                  boot_samp <- sample(plots, 20*ncopies, replace = T)
                  plots_samp <- sort(unique(boot_samp)) # the plots sampled
                  boot_table <- table(boot_samp) # a table of how often each plot was sampled
                  pn <- length(plots_samp) # the number of unique plots sampled
                  
                  # A loop to recreated the dataframe, by appending the dataset sampled plot at a time
                  for (p in 1:pn){
                    boot_df_i <- subset(Copy, plotBC == plots[p])
                    plot_n <- boot_table[p]
                    
                    k <- 1
                    while (plot_n >= k){
                      boot.df <- rbind(boot_df_i, boot.df)
                      k <- k+1
                    }
                  }
                  
                  n <- nrow(boot.df) # number of rows in bootstrapped dataframe
                  
                  #Step 2b. Simulate an intervention variable that is marginally independent of the simulated covariates
                  boot.df$X <- rbinom(n, 1, P.X)
                  boot.df$X1 <- 1
                  boot.df$X0 <- 0
                  
                  #Step 2c. Simulate each M as a function of its parents
                  boot.df$Mx <- get_m("X", binom = F)
                  boot.df$M1 <- get_m("X1", binom = F)
                  boot.df$M0 <- get_m("X0", binom = F)
                  
                  #Step 2d. simulate a Y for each type of effect as a function of its parents
                  
                  boot.df$Y_TE <- get_y("X", "Mx", binom = T)
                  boot.df$Y_PDE <- get_y("X", "M0", binom = T)
                  boot.df$Y_TIE <- get_y("X1", "Mx", binom = T)
                  boot.df$Y_TDE <- get_y("X", "M1", binom = T)
                  boot.df$Y_PIE <- get_y("X0", "Mx", binom = T)
                  
                  
                  #Step 3. Regress each potential health outcome against X
                  TE  <- coef(glm(Y_TE  ~ X, data = boot.df, family = "binomial"))[2]
                  PDE <- coef(glm(Y_PDE ~ X, data = boot.df, family = "binomial"))[2]
                  TIE <- coef(glm(Y_TIE ~ X, data = boot.df, family = "binomial"))[2]
                  TDE <- coef(glm(Y_TDE ~ X, data = boot.df, family = "binomial"))[2]
                  PIE <- coef(glm(Y_PIE ~ X, data = boot.df, family = "binomial"))[2]
                  
                  #Print summary
                  res <- c("TE" = TE, "PDE" = PDE, "TIE" = TIE, "TDE" = TDE, "PIE" = PIE)
                  
                  res
                }

stopCluster(cl)

## Summarize the results - median and 95% CI
MA_res <- matrixStats::colQuantiles(exp(res), probs = c(0.025, 0.5, 0.975), na.rm = T)
rownames(MA_res) <- c("TE", "PDE", "TIE", "TDE", "PIE")
MA_res

# percent mediated by burrows = TIE/TE
quantile(res[,3]/(res[,3] + res[,4])*100, probs = c(0.025, 0.5, 0.975), na.rm = T)
100-quantile(res[,3]/(res[,3] + res[,4])*100, probs = c(0.025, 0.5, 0.975), na.rm = T)


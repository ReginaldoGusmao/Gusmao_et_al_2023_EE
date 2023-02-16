# Gusmão & Gonçalves-souza 2021 ------------------------
  # Script 01_Route_full- Cap 1 
  # Main:
  #' SEM with comosite variables
  
  #--------------------------------------------------#
  
  # Starting R ####
# Clean memory
rm(list = ls())

# Directory

setwd("D:/OneDrive/Documentos/03_Doutorado/01_Artigo_trophic global/02_Capitulo_1_Diversity and body size/Git_Trophic_2.0/00_Data"
) 
# Data directory: https://1drv.ms/u/s!AnOMy2zrL6bDj6U5XrD11Y5G51ANTQ?e=Dk7ip7
# Password: Your nickname!

dir_pad <- getwd()

# Packages ####
# Linguagem
library(dplyr)
library(beepr)
# Analises
library(mgcv)
library(nlme)
library(vegan)
library(spdep)
library(adespatial)
library(piecewiseSEM)
# Graficos
library(ggplot2)
library(gridExtra)
library(corrplot)

#-----------------------------------------------------#
# Bird ####

## Import data ####
setwd(dir_pad)
setwd("anomaly")
anomaly_data <- read.csv("anomaly_bird_data.csv")

setwd(dir_pad)
setwd("Taxonomic")

bird <- read.csv("bird_cap1.csv",row.names = "X") 
bird <- cbind(bird,anomaly_data)

bird$median_BS_g <- log(bird$median_BS_g)


#-----------------------------------------------------#

## Create spacial Matrix ####
coords <- bird %>% dplyr::select("Longitude", "Latitude")

names(coords)
locations_sf <-
  st_as_sf(coords,
           coords = c("Longitude", "Latitude"),
           crs = 4326)
plot(locations_sf)

# Gerar lista W
coords_knn <-
  knearneigh(as.matrix(coords), k = 2, longlat = TRUE)
class(coords_knn)
coords_nb <- knn2nb(coords_knn, sym = TRUE)
coords_listw <- nb2listw(coords_nb, style = "W")
coords_listw

# dbMEM
MEM <-
  scores.listw(coords_listw, MEM.autocor = "positive")

candidates <-
  listw.candidates(
    coords,
    nb = c("gab", "mst", "dnear"),
    weights = c("binary", "flin")
  )

nbw <- length(candidates)
print(paste("Number of candidate SWMs generated:",nbw, sep = " "))

print(paste("New significance threshold value after p-value correction (Sidak correction)",round(1 - (1 - 0.05) ^ (1 / nbw), 4), sep = " "))

beep()

#-----------------------------------------------------#

## Richness model ####

### Create composite variable ####

# 1- observed model
ols_obs_ric <- lm(
  richness_bird ~  
    trophic_mean + median_BS_g + var_BS_g +
    bio_1  + bio_12 +
    anomaly_temp + anomaly_prec,
  data = bird
)

summary(ols_obs_ric)

# 2- Extract composite variables

summary(ols_obs_ric)$coefficients

### Traits composite
beta_trophic <- summary(ols_obs_ric)$coefficients[2,1]
  
beta_bs <- summary(ols_obs_ric)$coefficients[3,1]

beta_var <- summary(ols_obs_ric)$coefficients[4,1]

bird$composite_traits <- beta_bs * bird$median_BS_g +
  beta_var * bird$var_BS_g +
  beta_trophic * bird$trophic_mean

### Current Climate composite

beta_temp <- summary(ols_obs_ric)$coefficients[5,1]
  
beta_prec <- summary(ols_obs_ric)$coefficients[6,1]
  
bird$composite_cur_climate <- beta_temp * bird$bio_1 +
  beta_prec * bird$bio_12 

### Historical Climate composite

beta_an_temp <- summary(ols_obs_ric)$coefficients[7,1]
  
beta_an_prec <- summary(ols_obs_ric)$coefficients[8,1]
  
bird$composite_his_climate <- beta_an_temp * bird$anomaly_temp +
  beta_an_prec * bird$anomaly_prec

# correlação ####

cor_bird <- bird %>% 
  select("composite_cur_climate",
         "composite_his_climate",
         "bio_1",
         "bio_12",
         "anomaly_temp",
         "anomaly_prec")

cor_bird2 <- cor(cor_bird)
corrplot.mixed(cor_bird2)

cor_bird_2 <- bird %>% 
  select("composite_traits",
         "median_BS_g",
         "var_BS_g",
         "trophic_mean")
cor_bird2_2 <- cor(cor_bird_2)
corrplot.mixed(cor_bird2_2)

# 3- composite model
ols_comp_ric <- lm(
  richness_bird ~  
    composite_traits +
    composite_cur_climate +
    composite_his_climate,
  data = bird
)

summary(ols_comp_ric)

### Obtain MEM ####

res_ric <- ols_comp_ric$residuals

W_sel <-
  listw.select(
    res_ric,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_ric <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_ric, "MEMs_comp_bird_ric.csv", row.names = F)
MEMs_ric <- read.csv("MEMs_comp_bird_ric.csv")
setwd(dir_pad)

  # Add MEMs in df
bird_ric <- cbind(bird,MEMs_ric)

### OLS + MEMs ####
colnames(MEMs_ric)

mod_sp_comp_rich <- lm(
  richness_bird ~  
    composite_traits +
    composite_cur_climate +
    composite_his_climate+
    MEM5 + MEM4 + MEM8+
    MEM7 + MEM6 + MEM12 + MEM14,
  data = bird_ric
)

summary(mod_sp_comp_rich)

lm.morantest(mod_sp_comp_rich,coords_listw) # < .7
#-----------------------------------------------------#
## Composite traits ####
# 1- observed model

ols_obs_tra <- lm(
  composite_traits~  
    median_BS_g + var_BS_g +
    trophic_mean,
  data = bird
)

summary(ols_obs_tra)

### Obtain MEM ####

res_tra <- ols_obs_tra$residuals

W_sel <-
  listw.select(
    res_tra,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_tra <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_tra, "MEMs_comp_bird_tra_sim.csv", row.names = F)
MEMs_tra <- read.csv("MEMs_comp_bird_tra.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
bird_tra <- cbind(bird, MEMs_tra)

### OLS + MEMs ####
colnames(MEMs_tra)

mod_sp_comp_tra <-lm(
  composite_traits~  
    median_BS_g + var_BS_g +
    trophic_mean +
    MEM2,
  data = bird_tra
)

summary(mod_sp_comp_tra)

lm.morantest(mod_sp_comp_tra,coords_listw) # < .3


#-----------------------------------------------------#

## Trophic model ####

# 1- observed model
ols_obs_tro <- lm(
  trophic_mean ~  
    median_BS_g + var_BS_g +
    composite_cur_climate +
    composite_his_climate,
  data = bird
)

summary(ols_obs_tro)

### Obtain MEM ####

res_tro <- ols_obs_tro$residuals

W_sel <-
  listw.select(
    res_tro,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_tro <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_tro, "MEMs_comp_bird_tro.csv", row.names = F)
MEMs_tro <- read.csv("MEMs_comp_bird_tro.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
bird_tro <- cbind(bird, MEMs_tro)

### OLS + MEMs ####
colnames(MEMs_tro)

mod_sp_comp_tro <- lm(
  trophic_mean ~  
    composite_cur_climate +
    composite_his_climate +
    MEM1 + MEM6 + MEM2 +
    MEM7 + MEM8 + MEM12,
  data = bird_tro
)

summary(mod_sp_comp_tro)

lm.morantest(mod_sp_comp_tro,coords_listw) # < .5

#-----------------------------------------------------#

## Body size model ####

# 1- observed model
ols_obs_bs <- lm(
  median_BS_g ~ 
    composite_cur_climate +
    composite_his_climate,
  data = bird
)

summary(ols_obs_bs)

### Obtain MEM ####

res_bs <- ols_obs_bs$residuals

W_sel <-
  listw.select(
    res_bs,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_bs <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_bs, "MEMs_comp_bird_bs.csv", row.names = F)
MEMs_bs <- read.csv("MEMs_comp_bird_bs.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
bird_bs <- cbind(bird, MEMs_bs)

### OLS + MEMs ####
colnames(MEMs_bs)

mod_sp_comp_bs <- lm(
  median_BS_g ~  
    composite_cur_climate +
    composite_his_climate+
    MEM2 + MEM6  + MEM1 + MEM8 +
    MEM13+ MEM3 + MEM12 +MEM9,
  data = bird_bs
)

summary(mod_sp_comp_bs)

lm.morantest(mod_sp_comp_bs,coords_listw) # < .4

#-----------------------------------------------------#

## Variance Body size model ####

### Create composite variable ####

# 1- observed model
ols_obs_var_bs <- lm(
  var_BS_g ~ 
    composite_cur_climate +
    composite_his_climate,
  data = bird
)

summary(ols_obs_var_bs)


### Obtain MEM ####

res_var_bs <- ols_obs_var_bs$residuals

W_sel <-
  listw.select(
    res_var_bs,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_var_bs <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_var_bs, "MEMs_comp_bird_var_bs.csv", row.names = F)
MEMs_var_bs <- read.csv("MEMs_comp_bird_var_bs.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
bird_var_bs <- cbind(bird, MEMs_var_bs)

### OLS + MEMs ####
colnames(MEMs_var_bs)

mod_sp_comp_var_bs <- lm(
  var_BS_g ~
    composite_cur_climate +
    composite_his_climate +
    MEM2 + MEM1,
  data = bird_var_bs
)

summary(mod_sp_comp_var_bs)

lm.morantest(mod_sp_comp_var_bs,coords_listw) # < .8

#-----------------------------------------------------#
  
## Piecewise-SEM ####
mod_sem_com_bird <- piecewiseSEM::psem(mod_sp_comp_rich,
                                       mod_sp_comp_tra,
                                       mod_sp_comp_bs,
                                       mod_sp_comp_var_bs,
                                       mod_sp_comp_tro,
                                       median_BS_g %~~% var_BS_g,
                                       median_BS_g %~~%trophic_mean,
                                       var_BS_g %~~% trophic_mean)

dSep(mod_sem_com_bird)

beep()

sum_sem_com_bird<- summary(mod_sem_com_bird,
                           .progressBar = TRUE,
                           conserve = T)

beep()
sum_sem_com_bird


#-----------------------------------------------------#
# Mamm ####

## Import data ####
setwd(dir_pad)
setwd("anomaly")
anomaly_data <- read.csv("anomaly_mamm_data.csv")
setwd(dir_pad)

setwd(dir_pad)
setwd("Taxonomic")

mamm <- read.csv("mamm_cap1.csv",row.names = "X") 
setwd(dir_pad)

mamm <- cbind(mamm,anomaly_data)

mamm$median_BS_g <- log(mamm$median_BS_g)
#-----------------------------------------------------#

## Create spacial Matrix ####
coords <- mamm %>% dplyr::select("Longitude", "Latitude")

names(coords)
locations_sf <-
  st_as_sf(coords,
           coords = c("Longitude", "Latitude"),
           crs = 4326)
plot(locations_sf)

# Gerar lista W
coords_knn <-
  knearneigh(as.matrix(coords), k = 2, longlat = TRUE)
class(coords_knn)
coords_nb <- knn2nb(coords_knn, sym = TRUE)
coords_listw <- nb2listw(coords_nb, style = "W")
coords_listw

# dbMEM
MEM <-
  scores.listw(coords_listw, MEM.autocor = "positive")

candidates <-
  listw.candidates(
    coords,
    nb = c("gab", "mst", "dnear"),
    weights = c("binary", "flin")
  )

nbw <- length(candidates)
print(paste("Number of candidate SWMs generated:",nbw, sep = " "))

print(paste("New significance threshold value after p-value correction (Sidak correction)",round(1 - (1 - 0.05) ^ (1 / nbw), 4), sep = " "))

beep()

#-----------------------------------------------------#

## Richness model ####

### Create composite variable ####

# 1- observed model
ols_obs_ric <- lm(
  richness_mamm ~  
    trophic_mean + median_BS_g + var_BS_g +
    bio_1  + bio_12 +
    anomaly_temp + anomaly_prec,
  data = mamm
)

summary(ols_obs_ric)

# 2- Extract composite variables

summary(ols_obs_ric)$coefficients

### Traits composite
beta_trophic <- summary(ols_obs_ric)$coefficients[2,1]

beta_bs <- summary(ols_obs_ric)$coefficients[3,1]

beta_var <- summary(ols_obs_ric)$coefficients[4,1]

mamm$composite_traits <- beta_bs * mamm$median_BS_g +
  beta_var * mamm$var_BS_g +
  beta_trophic * mamm$trophic_mean

### Current Climate composite

beta_temp <- summary(ols_obs_ric)$coefficients[5,1]

beta_prec <- summary(ols_obs_ric)$coefficients[6,1]

mamm$composite_cur_climate <- beta_temp * mamm$bio_1 +
  beta_prec * mamm$bio_12 

### Historical Climate composite

beta_an_temp <- summary(ols_obs_ric)$coefficients[7,1]

beta_an_prec <- summary(ols_obs_ric)$coefficients[8,1]

mamm$composite_his_climate <- beta_an_temp * mamm$anomaly_temp +
  beta_an_prec * mamm$anomaly_prec 

# correlação ####

cor_mamm <- mamm %>% 
  select("composite_cur_climate",
         "composite_his_climate",
         "bio_1",
         "bio_12",
         "anomaly_temp",
         "anomaly_prec")

cor_mamm2 <- cor(cor_mamm)
corrplot.mixed(cor_mamm2)

cor_mamm_2 <- mamm %>% 
  select("composite_traits",
         "median_BS_g",
         "var_BS_g",
         "trophic_mean")
cor_mamm2_2 <- cor(cor_mamm_2)
corrplot.mixed(cor_mamm2_2)

# 3- composite model
ols_comp_ric <- lm(
  richness_mamm ~  
    composite_traits +
    composite_cur_climate +
    composite_his_climate,
  data = mamm
)

summary(ols_comp_ric)

### Obtain MEM ####

res_ric <- ols_comp_ric$residuals

W_sel <-
  listw.select(
    res_ric,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_ric <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_ric, "MEMs_comp_mamm_ric.csv", row.names = F)
MEMs_ric <- read.csv("MEMs_comp_mamm_ric.csv")
setwd(dir_pad)

# Add MEMs in df
mamm_ric <- cbind(mamm,MEMs_ric)

### OLS + MEMs ####
colnames(MEMs_ric)

mod_sp_comp_rich <- lm(
  richness_mamm ~  
    composite_traits +
    composite_cur_climate +
    composite_his_climate+
    MEM4 + MEM8 + MEM2 +
    MEM7 + MEM5 + MEM1 + 
    MEM3 + MEM6 + MEM12 +
    MEM11,
  data = mamm_ric
)

summary(mod_sp_comp_rich)

lm.morantest(mod_sp_comp_rich,coords_listw) # < .7
#-----------------------------------------------------#

## Composite traits ####
# 1- observed model

ols_obs_tra <- lm(
  composite_traits~  
    median_BS_g + var_BS_g +
    trophic_mean,
  data = mamm
)

summary(ols_obs_tra)

### Obtain MEM ####

res_tra <- ols_obs_tra$residuals

W_sel <-
  listw.select(
    res_tra,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_tra <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_tra, "MEMs_comp_mamm_tra_sim.csv", row.names = F)
MEMs_tra <- read.csv("MEMs_comp_mamm_tra_sim.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
mamm_tra <- cbind(mamm, MEMs_tra)

### OLS + MEMs ####
colnames(MEMs_tra)

mod_sp_comp_tra <-lm(
  composite_traits~  
    median_BS_g + var_BS_g +
    trophic_mean +
    MEM51 + MEM53 +
    MEM54 + MEM57,
  data = mamm_tra
)

summary(mod_sp_comp_tra)

lm.morantest(mod_sp_comp_tra,coords_listw) # < .7

#-----------------------------------------------------#

## Trophic model ####

# 1- observed model
ols_obs_tro <- lm(
  trophic_mean ~  
    composite_cur_climate +
    composite_his_climate,
  data = mamm
)

summary(ols_obs_tro)

### Obtain MEM ####

res_tro <- ols_obs_tro$residuals

W_sel <-
  listw.select(
    res_tro,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_tro <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_tro, "MEMs_comp_mamm_tro.csv", row.names = F)
MEMs_tro <- read.csv("MEMs_comp_mamm_tro.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
mamm_tro <- cbind(mamm, MEMs_tro)

### OLS + MEMs ####
colnames(MEMs_tro)

mod_sp_comp_tro <- lm(
  trophic_mean ~  
    composite_cur_climate +
    composite_his_climate +
    MEM6 + MEM4 + MEM3 +
    MEM7 + MEM15 + MEM21,
  data = mamm_tro
)

summary(mod_sp_comp_tro)

lm.morantest(mod_sp_comp_tro,coords_listw) # < .7

#-----------------------------------------------------#

## Body size model ####

# 1- observed model
ols_obs_bs <- lm(
  median_BS_g ~ 
    composite_cur_climate +
    composite_his_climate,
  data = mamm
)

summary(ols_obs_bs)

### Obtain MEM ####

res_bs <- ols_obs_bs$residuals

W_sel <-
  listw.select(
    res_bs,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_bs <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_bs, "MEMs_comp_mamm_bs.csv", row.names = F)
MEMs_bs <- read.csv("MEMs_comp_mamm_bs.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
mamm_bs <- cbind(mamm, MEMs_bs)

### OLS + MEMs ####
colnames(MEMs_bs)

mod_sp_comp_bs <- lm(
  median_BS_g ~  
    composite_cur_climate +
    composite_his_climate+
    MEM4 + MEM1  + MEM3 + MEM6 +
    MEM7+ MEM8 + MEM10 + MEM18 +
    MEM13 + MEM14,
  data = mamm_bs
)

summary(mod_sp_comp_bs)

lm.morantest(mod_sp_comp_bs,coords_listw) # < .7

#-----------------------------------------------------#
  
 ## Variance Body size model ####

### Create composite variable ####

# 1- observed model
ols_obs_var_bs <- lm(
  var_BS_g ~ 
    composite_cur_climate +
    composite_his_climate,
  data = mamm
)

summary(ols_obs_var_bs)

### Obtain MEM ####

res_var_bs <- ols_obs_var_bs$residuals

W_sel <-
  listw.select(
    res_var_bs,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_var_bs <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_var_bs, "MEMs_comp_mamm_var_bs.csv", row.names = F)
MEMs_var_bs <- read.csv("MEMs_comp_mamm_var_bs.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
mamm_var_bs <- cbind(mamm, MEMs_var_bs)

### OLS + MEMs ####
colnames(MEMs_var_bs)

mod_sp_comp_var_bs <- lm(
  var_BS_g ~
    composite_cur_climate +
    composite_his_climate +
    MEM1 + MEM4 + MEM5 + 
    MEM8 + MEM7 + MEM3 + MEM6,
  data = mamm_var_bs
)

summary(mod_sp_comp_var_bs)

lm.morantest(mod_sp_comp_var_bs,coords_listw) # < .7

#-----------------------------------------------------#
## Piecewise-SEM ####

mod_sem_com_mamm <- piecewiseSEM::psem(mod_sp_comp_rich,
                                       mod_sp_comp_tra,
                                       mod_sp_comp_bs,
                                       mod_sp_comp_var_bs,
                                       mod_sp_comp_tro,
                                       median_BS_g %~~% var_BS_g,
                                       median_BS_g %~~%trophic_mean,
                                       var_BS_g %~~% trophic_mean)

dSep(mod_sem_com_mamm)

beep()

sum_sem_com_mamm<- summary(mod_sem_com_mamm,
                           .progressBar = TRUE,
                           conserve = T)

beep()
sum_sem_com_mamm

#-----------------------------------------------------#
# Amph ####

## Import data ####
setwd(dir_pad)
setwd("anomaly")
anomaly_data <- read.csv("anomaly_amph_data.csv")
setwd(dir_pad)

setwd(dir_pad)
setwd("Taxonomic")

amph <- read.csv("amph_cap1.csv",row.names = "X") 
setwd(dir_pad)

amph <- cbind(amph,anomaly_data)

amph$median_BS_g <- log(amph$median_BS_g)



#-----------------------------------------------------#

## Create spacial Matrix ####
coords <- amph %>% dplyr::select("Longitude", "Latitude")

names(coords)
locations_sf <-
  st_as_sf(coords,
           coords = c("Longitude", "Latitude"),
           crs = 4326)
plot(locations_sf)

# Gerar lista W
coords_knn <-
  knearneigh(as.matrix(coords), k = 2, longlat = TRUE)
class(coords_knn)
coords_nb <- knn2nb(coords_knn, sym = TRUE)
coords_listw <- nb2listw(coords_nb, style = "W")
coords_listw

# dbMEM
MEM <-
  scores.listw(coords_listw, MEM.autocor = "positive")

candidates <-
  listw.candidates(
    coords,
    nb = c("gab", "mst", "dnear"),
    weights = c("binary", "flin")
  )

nbw <- length(candidates)
print(paste("Number of candidate SWMs generated:",nbw, sep = " "))

print(paste("New significance threshold value after p-value correction (Sidak correction)",round(1 - (1 - 0.05) ^ (1 / nbw), 4), sep = " "))

beep()

#-----------------------------------------------------#

## Richness model ####

### Create composite variable ####

# 1- observed model
ols_obs_ric <- lm(
  richness_amph ~  
    trophic_mean + median_BS_g + var_BS_g +
    bio_1  + bio_12 +
    anomaly_prec, #+
    #anomaly_temp,
  data = amph
)

summary(ols_obs_ric)

# 2- Extract composite variables

summary(ols_obs_ric)$coefficients

### Traits composite
beta_trophic <- summary(ols_obs_ric)$coefficients[2,1]

beta_bs <- summary(ols_obs_ric)$coefficients[3,1]

beta_var <- summary(ols_obs_ric)$coefficients[4,1]

amph$composite_traits <- beta_bs * amph$median_BS_g +
  beta_var * amph$var_BS_g +
  beta_trophic * amph$trophic_mean

### Current Climate composite

beta_temp <- summary(ols_obs_ric)$coefficients[5,1]

beta_prec <- summary(ols_obs_ric)$coefficients[6,1]

amph$composite_cur_climate <- beta_temp * amph$bio_1 +
  beta_prec * amph$bio_12 

##### Historical Climate composite

beta_an_temp <- summary(ols_obs_ric)$coefficients[7,1]

#beta_an_prec <- summary(ols_obs_ric)$coefficients[8,1]

amph$composite_his_climate <- beta_an_temp * amph$anomaly_temp +
  beta_an_prec * amph$anomaly_prec 

# correlação ####

cor_amph <- amph %>% 
  select("composite_cur_climate",
         "composite_his_climate",
         "bio_1",
         "bio_12",
         "anomaly_temp",
         "anomaly_prec")

cor_amph2 <- cor(cor_amph)
corrplot.mixed(cor_amph2)

cor_amph_2 <- amph %>% 
  select("composite_traits",
         "median_BS_g",
         "var_BS_g",
         "trophic_mean")
cor_amph2_2 <- cor(cor_amph_2)
corrplot.mixed(cor_amph2_2)

# 3- composite model
ols_comp_ric <- lm(
  richness_amph~  
    composite_traits +
    composite_cur_climate +
    anomaly_prec,
  data = amph
)

summary(ols_comp_ric)

### Obtain MEM ####

res_ric <- ols_comp_ric$residuals

W_sel <-
  listw.select(
    res_ric,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_ric <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_ric, "MEMs_comp_amph_ric.csv", row.names = F)
MEMs_ric <- read.csv("MEMs_comp_amph_ric.csv")
setwd(dir_pad)

# Add MEMs in df
amph_ric <- cbind(amph,MEMs_ric)

### OLS + MEMs ####
colnames(MEMs_ric)

mod_sp_comp_rich <- lm(
  richness_amph ~  
    composite_traits +
    composite_cur_climate +
    anomaly_prec+
MEM1 + MEM3 + MEM2+
  MEM5 + MEM4,
  data = amph_ric
)

summary(mod_sp_comp_rich)

lm.morantest(mod_sp_comp_rich,coords_listw) # < .7
#-----------------------------------------------------#

## Composite traits ####
# 1- observed model

ols_obs_tra <- lm(
  composite_traits~  
    median_BS_g + var_BS_g +
    trophic_mean,
  data = amph
)

summary(ols_obs_tra)

### Obtain MEM ####

res_tra <- ols_obs_tra$residuals

W_sel <-
  listw.select(
    res_tra,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_tra <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_tra, "MEMs_comp_amph_tra_sim.csv", row.names = F)
MEMs_tra <- read.csv("MEMs_comp_amph_tra_sim.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
amph_tra <- cbind(amph, MEMs_tra)

### OLS + MEMs ####
colnames(MEMs_tra)

mod_sp_comp_tra <-lm(
  composite_traits~  
    median_BS_g + var_BS_g +
    trophic_mean +
 MEM401,
  data = amph_tra
)

summary(mod_sp_comp_tra)

lm.morantest(mod_sp_comp_tra,coords_listw) # < -.003

#-----------------------------------------------------#

## Trophic model ####

# 1- observed model
ols_obs_tro <- lm(
  trophic_mean ~  
    composite_cur_climate +
    anomaly_prec,
  data = amph
)

summary(ols_obs_tro)

### Obtain MEM ####

res_tro <- ols_obs_tro$residuals

W_sel <-
  listw.select(
    res_tro,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_tro <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_tro, "MEMs_comp_amph_tro.csv", row.names = F)
MEMs_tro <- read.csv("MEMs_comp_amph_tro.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
amph_tro <- cbind(amph, MEMs_tro)

### OLS + MEMs ####
colnames(MEMs_tro)

mod_sp_comp_tro <- lm(
  trophic_mean ~  
    composite_cur_climate +
    anomaly_prec +
MEM2 +MEM5 + MEM3 + MEM4 +
  MEM7 + MEM10 + MEM14 +MEM1 +
  MEM13,
  data = amph_tro
)

summary(mod_sp_comp_tro)

lm.morantest(mod_sp_comp_tro,coords_listw) # < .7

#-----------------------------------------------------#

## Body size model ####

# 1- observed model
ols_obs_bs <- lm(
  median_BS_g ~ 
    composite_cur_climate +
    anomaly_prec,
  data = amph
)

summary(ols_obs_bs)

### Obtain MEM ####

res_bs <- ols_obs_bs$residuals

W_sel <-
  listw.select(
    res_bs,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_bs <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_bs, "MEMs_comp_amph_bs.csv", row.names = F)
MEMs_bs <- read.csv("MEMs_comp_amph_bs.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
amph_bs <- cbind(amph, MEMs_bs)

### OLS + MEMs ####
colnames(MEMs_bs)

mod_sp_comp_bs <- lm(
  median_BS_g ~  
    composite_cur_climate +
    anomaly_prec+
  MEM5 + MEM4 + MEM8 + MEM3 +
    MEM2 + MEM12+ MEM7 +MEM6,
  data = amph_bs
)

summary(mod_sp_comp_bs)

lm.morantest(mod_sp_comp_bs,coords_listw) # < .6

#-----------------------------------------------------#

## Variance Body size model ####

### Create composite variable ####

# 1- observed model
ols_obs_var_bs <- lm(
  var_BS_g ~ 
    composite_cur_climate +
    anomaly_prec,
  data = amph
)

summary(ols_obs_var_bs)

### Obtain MEM ####

res_var_bs <- ols_obs_var_bs$residuals

W_sel <-
  listw.select(
    res_var_bs,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_var_bs <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_var_bs, "MEMs_comp_amph_var_bs.csv", row.names = F)
MEMs_var_bs <- read.csv("MEMs_comp_amph_var_bs.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
amph_var_bs <- cbind(amph, MEMs_var_bs)

### OLS + MEMs ####
colnames(MEMs_var_bs)

mod_sp_comp_var_bs <- lm(
  var_BS_g ~
    composite_cur_climate +
    anomaly_prec +
    MEM7,
  data = amph_var_bs
)

summary(mod_sp_comp_var_bs)

lm.morantest(mod_sp_comp_var_bs,coords_listw) # < .4

#-----------------------------------------------------#
# Piecewise-SEM ####

mod_sem_com_amph <- piecewiseSEM::psem(mod_sp_comp_rich,
                                       mod_sp_comp_tra,
                                       mod_sp_comp_bs,
                                       mod_sp_comp_var_bs,
                                       mod_sp_comp_tro,
                                       median_BS_g %~~% var_BS_g,
                                       median_BS_g %~~%trophic_mean,
                                       var_BS_g %~~% trophic_mean)

dSep(mod_sem_com_amph)

beep()

sum_sem_com_amph<- summary(mod_sem_com_amph,
                           .progressBar = TRUE,
                           conserve = T)

beep()
sum_sem_com_amph


#-----------------------------------------------------#

# Lizard ####

## Import data ####
setwd(dir_pad)
setwd("anomaly")
anomaly_data <- read.csv("anomaly_lizard_data.csv")
setwd(dir_pad)

setwd(dir_pad)
setwd("Taxonomic")

lizard <- read.csv("lizard_cap1.csv",row.names = "X") 
setwd(dir_pad)

lizard <- cbind(lizard,anomaly_data)

lizard$median_BS_g <- log(lizard$median_BS_g)
#-----------------------------------------------------#

## Create spacial Matrix ####
coords <- lizard %>% dplyr::select("Longitude", "Latitude")

names(coords)
locations_sf <-
  st_as_sf(coords,
           coords = c("Longitude", "Latitude"),
           crs = 4326)
plot(locations_sf)

# Gerar lista W
coords_knn <-
  knearneigh(as.matrix(coords), k = 2, longlat = TRUE)
class(coords_knn)
coords_nb <- knn2nb(coords_knn, sym = TRUE)
coords_listw <- nb2listw(coords_nb, style = "W")
coords_listw

# dbMEM
MEM <-
  scores.listw(coords_listw, MEM.autocor = "positive")

candidates <-
  listw.candidates(
    coords,
    nb = c("gab", "mst", "dnear"),
    weights = c("binary", "flin")
  )

nbw <- length(candidates)
print(paste("Number of candidate SWMs generated:",nbw, sep = " "))

print(paste("New significance threshold value after p-value correction (Sidak correction)",round(1 - (1 - 0.05) ^ (1 / nbw), 4), sep = " "))

beep()

#-----------------------------------------------------#

## Richness model ####

### Create composite variable ####

# 1- observed model
ols_obs_ric <- lm(
  richness_lizard ~  
    trophic_mean + median_BS_g + var_BS_g +
    bio_1  + bio_12 +
    #anomaly_temp + 
    anomaly_prec,
  data = lizard
)

summary(ols_obs_ric)

# 2- Extract composite variables

summary(ols_obs_ric)$coefficients

### Traits composite
beta_trophic <- summary(ols_obs_ric)$coefficients[2,1]

beta_bs <- summary(ols_obs_ric)$coefficients[3,1]

beta_var <- summary(ols_obs_ric)$coefficients[4,1]

lizard$composite_traits <- beta_bs * lizard$median_BS_g +
  beta_var * lizard$var_BS_g +
  beta_trophic * lizard$trophic_mean

### Current Climate composite

beta_temp <- summary(ols_obs_ric)$coefficients[5,1]

beta_prec <- summary(ols_obs_ric)$coefficients[6,1]

lizard$composite_cur_climate <- beta_temp * lizard$bio_1 +
  beta_prec * lizard$bio_12 

##### Historical Climate composite

beta_an_temp <- summary(ols_obs_ric)$coefficients[7,1]

#beta_an_prec <- summary(ols_obs_ric)$coefficients[8,1]

lizard$composite_his_climate <- beta_an_temp * lizard$anomaly_temp +
  beta_an_prec * lizard$anomaly_prec

# correlação ####

cor_lizard <- lizard %>% 
  select("composite_cur_climate",
         "composite_his_climate",
         "bio_1",
         "bio_12",
         "anomaly_temp",
         "anomaly_prec")

cor_lizard2 <- cor(cor_lizard)
corrplot.mixed(cor_lizard2)

cor_lizard_2 <- lizard %>% 
  select("composite_traits",
         "median_BS_g",
         "var_BS_g",
         "trophic_mean")
cor_lizard2_2 <- cor(cor_lizard_2)
corrplot.mixed(cor_lizard2_2)

# 3- composite model
ols_comp_ric <- lm(
  richness_lizard~  
    composite_traits +
    composite_cur_climate +
    anomaly_prec,
  data = lizard
)

summary(ols_comp_ric)

### Obtain MEM ####

res_ric <- ols_comp_ric$residuals

W_sel <-
  listw.select(
    res_ric,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_ric <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_ric, "MEMs_comp_lizard_ric.csv", row.names = F)
MEMs_ric <- read.csv("MEMs_comp_lizard_ric.csv")
setwd(dir_pad)

# Add MEMs in df
lizard_ric <- cbind(lizard,MEMs_ric)

### OLS + MEMs ####
colnames(MEMs_ric)

mod_sp_comp_rich <- lm(
  richness_lizard ~  
    composite_traits +
    composite_cur_climate +
    anomaly_prec+
MEM5 + MEM1 + MEM3 +MEM6 +MEM8 + 
  MEM12 + MEM11 + MEM7 +MEM10 +MEM16 +
  MEM18,
  data = lizard_ric
)

summary(mod_sp_comp_rich)

lm.morantest(mod_sp_comp_rich,coords_listw) # < .7
#-----------------------------------------------------#

## Composite traits ####
# 1- observed model

ols_obs_tra <- lm(
  composite_traits~  
    median_BS_g + var_BS_g +
    trophic_mean,
  data = lizard
)

summary(ols_obs_tra)

### Obtain MEM ####

res_tra <- ols_obs_tra$residuals

W_sel <-
  listw.select(
    res_tra,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_tra <- W_sel$best$MEM.select

###############-

# Add MEMs in df
lizard_tra <- lizard

### OLS + MEMs ####
colnames(MEMs_tra)

mod_sp_comp_tra <-lm(
  composite_traits~  
    median_BS_g + var_BS_g +
    trophic_mean,
  data = lizard_tra
)

summary(mod_sp_comp_tra)

lm.morantest(mod_sp_comp_tra,coords_listw) # < -.010

#-----------------------------------------------------#
## Trophic model ####

# 1- observed model
ols_obs_tro <- lm(
  trophic_mean ~  
    composite_cur_climate +
    anomaly_prec,
  data = lizard
)

summary(ols_obs_tro)

### Obtain MEM ####

res_tro <- ols_obs_tro$residuals

W_sel <-
  listw.select(
    res_tro,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_tro <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_tro, "MEMs_comp_lizard_tro.csv", row.names = F)
MEMs_tro <- read.csv("MEMs_comp_lizard_tro.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
lizard_tro <- cbind(lizard, MEMs_tro)

### OLS + MEMs ####
colnames(MEMs_tro)

mod_sp_comp_tro <- lm(
  trophic_mean ~  
    composite_cur_climate +
    anomaly_prec +
MEM4 + MEM6 + MEM3 + MEM2 + MEM9 +
  MEM7 + MEM12 + MEM10,
  data = lizard_tro
)

summary(mod_sp_comp_tro)

lm.morantest(mod_sp_comp_tro,coords_listw) # < .7

#-----------------------------------------------------#

## Body size model ####

# 1- observed model
ols_obs_bs <- lm(
  median_BS_g ~ 
    composite_cur_climate +
    anomaly_prec,
  data = lizard)

summary(ols_obs_bs)

### Obtain MEM ####

res_bs <- ols_obs_bs$residuals

W_sel <-
  listw.select(
    res_bs,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_bs <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_bs, "MEMs_comp_lizard_bs.csv", row.names = F)
MEMs_bs <- read.csv("MEMs_comp_lizard_bs.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
lizard_bs <- cbind(lizard, MEMs_bs)

### OLS + MEMs ####
colnames(MEMs_bs)

mod_sp_comp_bs <- lm(
  median_BS_g ~  
    composite_cur_climate +
    anomaly_prec+
 MEM4 + MEM1 +MEM2 + MEM6 + 
   MEM3 +MEM9 + MEM11,
  data = lizard_bs
)

summary(mod_sp_comp_bs)

lm.morantest(mod_sp_comp_bs,coords_listw) # < .5

#-----------------------------------------------------#

## Variance Body size model ####

### Create composite variable ####

# 1- observed model
ols_obs_var_bs <- lm(
  var_BS_g ~ 
    composite_cur_climate +
    anomaly_prec,
  data = lizard
)

summary(ols_obs_var_bs)

### Obtain MEM ####

res_var_bs <- ols_obs_var_bs$residuals

W_sel <-
  listw.select(
    res_var_bs,
    candidates,
    MEM.autocor = "positive",
    p.adjust = TRUE,
    method = "MIR",
    verbose = TRUE
  )
beep()

W_sel$best.id

W_sel$best$MEM.select
W_sel$candidates$N.var[W_sel$best.id]
W_sel$candidates$Pvalue[W_sel$best.id]

#Final spatial matrix
MEMs_var_bs <- W_sel$best$MEM.select

#####  Save MEMS ##### 
setwd(dir_pad)
setwd("MEMS/01_Taxonomic")
#write.csv(MEMs_var_bs, "MEMs_comp_lizard_var_bs.csv", row.names = F)
MEMs_var_bs <- read.csv("MEMs_comp_lizard_var_bs.csv")
setwd(dir_pad)

###############-

# Add MEMs in df
lizard_var_bs <- cbind(lizard, MEMs_var_bs)

### OLS + MEMs ####
colnames(MEMs_var_bs)

mod_sp_comp_var_bs <- lm(
  var_BS_g ~
    composite_cur_climate +
    anomaly_prec +
    MEM2 + MEM3 + MEM5 + MEM11 +MEM13 + MEM4 + MEM10 +
    MEM17 + MEM14 + MEM6 + MEM8 + MEM18,
  data = lizard_var_bs
)

summary(mod_sp_comp_var_bs)

lm.morantest(mod_sp_comp_var_bs,coords_listw) # < .4

#-----------------------------------------------------#
# Piecewise-SEM ####

mod_sem_com_lizard <- piecewiseSEM::psem(mod_sp_comp_rich,
                                       mod_sp_comp_tra,
                                       mod_sp_comp_bs,
                                       mod_sp_comp_var_bs,
                                       mod_sp_comp_tro,
                                       median_BS_g %~~% var_BS_g,
                                       median_BS_g %~~%trophic_mean,
                                       var_BS_g %~~% trophic_mean)

dSep(mod_sem_com_lizard)

beep()

sum_sem_com_lizard<- summary(mod_sem_com_lizard,
                           .progressBar = TRUE,
                           conserve = T)

beep()
sum_sem_com_lizard









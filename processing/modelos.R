# Librerías básicas
library(tidyverse)
library(sjlabelled)

# Librerias análisis de secuencia
library(TraMineR)
library(TraMineRextras)
library(WeightedCluster)
library(cluster)

# Librería ELSOC
library(elsoc)

# Librerias graficos
library(ggpubr)
library(shiny)
library(plotly)
library(ggstance)

# Libreria analisis estadístico
library(nnet)
library(modelsummary)
library(lcmm)
library(rcompanion)
library(statar)
library(marginaleffects)

set.seed(123)


# Cargar datos ELSOC:
# elsoc::load_elsoc(data = 'wide')
load(file.path('..', 'inputs', 'ELSOC_Wide_2016_2021_v1.00_R.RData'))

# Crear variables de sintomatología depresiva por ola:
elsoc_salud <- elsoc_wide_2016_2021 %>% 
  
  # Se conservan observaciones presentes en las olas 1 a 4 (muestra original)
  filter(tipo_atricion %in% 1:2) %>% 
  
  # Crear variables de PHQ9
  purrr::map_at(.at = vars(starts_with('s11_0')), 
                .f = function(s) car::recode(s, "1 = 0; 2 = 1; 3 = 2; c(4, 5) = 3; c(-888, -999) = NA")) %>%
  as.data.frame() %>%  
  mutate(phq9_w01 = (s11_01_w01 + s11_02_w01 + s11_03_w01 + s11_04_w01 + s11_05_w01 + s11_06_w01 + s11_07_w01 + s11_08_w01 + s11_09_w01),
         phq9_w02 = (s11_01_w02 + s11_02_w02 + s11_03_w02 + s11_04_w02 + s11_05_w02 + s11_06_w02 + s11_07_w02 + s11_08_w02 + s11_09_w02),
         phq9_w03 = (s11_01_w03 + s11_02_w03 + s11_03_w03 + s11_04_w03 + s11_05_w03 + s11_06_w03 + s11_07_w03 + s11_08_w03 + s11_09_w03),
         phq9_w04 = (s11_01_w04 + s11_02_w04 + s11_03_w04 + s11_04_w04 + s11_05_w04 + s11_06_w04 + s11_07_w04 + s11_08_w04 + s11_09_w04)
  ) %>% 
  
  # Quitar NAs
  filter(!is.na(phq9_w01), !is.na(phq9_w02), !is.na(phq9_w03), !is.na(phq9_w04)) %>% 
  
  # Crear indicador de depresión en 4 y en 2 categorías:
  mutate(depr4_w01 = car::recode(phq9_w01, "0:4 = 1; 5:9 = 2; 10:14 = 3; 15:27 = 4"), 
         depr4_w02 = car::recode(phq9_w02, "0:4 = 1; 5:9 = 2; 10:14 = 3; 15:27 = 4"), 
         depr4_w03 = car::recode(phq9_w03, "0:4 = 1; 5:9 = 2; 10:14 = 3; 15:27 = 4"), 
         depr4_w04 = car::recode(phq9_w04, "0:4 = 1; 5:9 = 2; 10:14 = 3; 15:27 = 4"), 
         depr2_w01 = car::recode(phq9_w01, "0:9 = 0; 10:27 = 1"), 
         depr2_w02 = car::recode(phq9_w02, "0:9 = 0; 10:27 = 1"), 
         depr2_w03 = car::recode(phq9_w03, "0:9 = 0; 10:27 = 1"),  
         depr2_w04 = car::recode(phq9_w04, "0:9 = 0; 10:27 = 1"), 
         phq_sum = depr4_w01 + depr4_w02 + depr4_w03 + depr4_w04
  )


# Cargar datos ELSOC LONG:
# elsoc::load_elsoc(data = 'long')
load(file.path('..', 'inputs', 'ELSOC_Long_2016_2021_v1.00_R.RData'))

# Crear variables de sintomatología depresiva por ola LONG:
elsoc_salud_long <- elsoc_long_2016_2021 %>% 
  
  # Se conservan observaciones presentes en las 5 olas
  filter(tipo_atricion %in% 1:2, ola != 5) %>% 
  
  # Crear variables de PHQ9
  purrr::map_at(.at = vars(starts_with('s11_0')), 
                .f = function(s) car::recode(s, "1 = 0; 2 = 1; 3 = 2; c(4, 5) = 3; c(-888, -999) = NA")) %>%
  as.data.frame() %>%  
  mutate(phq9 = (s11_01 + s11_02 + s11_03 + s11_04 + s11_05 + s11_06 + s11_07 + s11_08 + s11_09)) %>% 
  
  # Quitar NAs
  filter(!is.na(phq9)) %>% 
  
  # Crear indicador de depresión en 4 y en 2 categorías:
  mutate(depr4 = car::recode(phq9, "0:4 = 1; 5:9 = 2; 10:14 = 3; 15:27 = 4"), 
         depr2 = car::recode(phq9, "0:9 = 0; 10:27 = 1"))



# Se usa el modelo de 1 clase para fijar los valores iniciales a iterar:
  

lcmm1_lin <- hlme(phq9 ~ ola, subject = "idencuesta", ng = 1,  
                  data = elsoc_salud_long)

lcmm3_lin <- gridsearch(rep = 5, maxiter = 5, minit = lcmm1_lin, 
                        hlme(phq9 ~ ola, mixture = ~ ola, subject = "idencuesta", ng = 3, 
                             data = elsoc_salud_long))

lcmm4_lin <- gridsearch(rep = 5, maxiter = 5, minit = lcmm1_lin, 
                        hlme(phq9 ~ ola, mixture = ~ ola,  subject = "idencuesta", ng = 4, 
                             data = elsoc_salud_long))

lcmm5_lin <- gridsearch(rep = 5, maxiter = 5, minit = lcmm1_lin, 
                        hlme(phq9 ~ ola, mixture = ~ ola,  subject = "idencuesta", ng = 5, 
                             data = elsoc_salud_long))

save(lcmm1_lin, lcmm3_lin, lcmm4_lin, lcmm5_lin, 
     file = file.path('..', 'inputs', 'modelos_lcmm.RData'))


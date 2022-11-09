# Librerías básicas
library(tidyverse)
library(sjlabelled)

# Librería ELSOC
library(elsoc)

# Libreria analisis estadístico
library(lcmm)

set.seed(123)

# Cargar datos ELSOC:
load(file.path('..', 'inputs', 'datos_elsoc_preparados.RData'))

# Se usa el modelo de 1 clase para fijar los valores iniciales a iterar:

lcmm1_lin <- hlme(phq9 ~ ola, subject = "idencuesta", ng = 1,  
                  data = elsoc_salud_long)

lcmm3_lin <- gridsearch(rep = 5, maxiter = 5, minit = lcmm1_lin, 
                        hlme(phq9 ~ ola, mixture = ~ ola, 
                             subject = "idencuesta", ng = 3, 
                             data = elsoc_salud_long))

lcmm4_lin <- gridsearch(rep = 5, maxiter = 5, minit = lcmm1_lin, 
                        hlme(phq9 ~ ola, mixture = ~ ola,  
                             subject = "idencuesta", ng = 4, 
                             data = elsoc_salud_long))

lcmm5_lin <- gridsearch(rep = 5, maxiter = 5, minit = lcmm1_lin, 
                        hlme(phq9 ~ ola, mixture = ~ ola,  
                             subject = "idencuesta", ng = 5, 
                             data = elsoc_salud_long))


elsoc_lcmm <- elsoc_salud_wide %>%
  select(idencuesta) %>% 
  left_join(lcmm4_lin$pprob %>% select(idencuesta, class), 
            by = 'idencuesta') %>% 
  rename(class4 = class) %>% 
  left_join(lcmm5_lin$pprob %>% select(idencuesta, class),
            by = 'idencuesta') %>%
  rename(class5 = class) %>% 
  mutate(trayectoria_gmm.4lin = factor(class4,
                                       levels = c(2, 1, 3, 4),
                                       labels = c('Trayectoria asintomatica',
                                                  'Trayectoria sintomas moderados',
                                                  'Trayectoria de sintomas ascendente',
                                                  'Trayectoria sintomas severos persistentes')),
         trayectoria_gmm.5lin = factor(class5, levels = c(4, 5, 1, 2, 3),
                                       labels = c('Trayectoria baja-media',
                                                  'Trayectoria media-moderada',
                                                  'Trayectoria descendente',
                                                  'Trayectoria ascendente',
                                                  'Trayectoria alta-moderada')))
         
save(elsoc_lcmm, lcmm1_lin, lcmm3_lin, lcmm4_lin, lcmm5_lin, 
     file = file.path('..', 'inputs', 'modelos_lcmm.RData'))


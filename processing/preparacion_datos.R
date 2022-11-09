# Librerías básicas
library(tidyverse)
library(sjlabelled)

# Librería ELSOC
library(elsoc)

###################################
### Datos ELSOC en formato WIDE ###
###################################

elsoc::load_elsoc(data = 'wide')


tramos_m30 <- "1 = 220; 2 = 250; 3 = 305; 4 = 355; 5 = 400; 6 = 445; 7 = 490; 8 = 535; 9 = 585; 10 = 640; 11 = 700; 12 = 765; 13 = 845; 14 = 935; 15 = 1040; 16 = 1180; 17 = 1375; 18 = 1670; 19 = 2275; 20 = 2700; else = NA"

elsoc_salud_wide <- elsoc_wide_2016_2021 %>% 
  
  # Se conservan observaciones completas presentes en las olas 1 a 4 (muestra original)
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
  ) %>% 
  
  # Variables adicionales:

  mutate(
    
    # Variables de ingreso:
    m30_w01 = 1000*as.numeric(car::recode(m30_w01, tramos_m30)),
    m29_w01_imp = ifelse(is_nsnr(m29_w01), m30_w01, m29_w01),
    ypc_w01 = m29_w01_imp / nhogar1_w01,
    
    m30_w02 = 1000*as.numeric(car::recode(m30_w02, tramos_m30)),
    m29_w02_imp = ifelse(is_nsnr(m29_w02), m30_w02, m29_w02),
    ypc_w02 = m29_w02_imp / nhogar1_w01,    
    
    m30_w03 = 1000*as.numeric(car::recode(m30_w03, tramos_m30)),
    m29_w03_imp = ifelse(is_nsnr(m29_w03), m30_w03, m29_w03),
    ypc_w03 = m29_w03_imp / ifelse(is_nsnr(m54_w03), NA, m54_w03),
    
    m30_w04 = 1000*as.numeric(car::recode(m30_w04, tramos_m30)),
    m29_w04_imp = ifelse(is_nsnr(m29_w04), m30_w04, m29_w04),
    ypc_w04 = m29_w04_imp / ifelse(is_nsnr(m54_w04), NA, m54_w04),
    
    log_ing = log((ypc_w01 + ypc_w02 + ypc_w03 + ypc_w04)/5),
    
    quintil_w01 = statar::xtile(ypc_w01, n = 5, wt = ponderador02_w01),
    quintil_w01 = factor(quintil_w01,
                         levels = 1:5,
                         labels = glue::glue('Quintil {1:5}')),
    quintil_w02 = statar::xtile(ypc_w02, n = 5, wt = ponderador02_w02),
    quintil_w02 = factor(quintil_w02,
                         levels = 1:5,
                         labels = glue::glue('Quintil {1:5}')),
    quintil_w03 = statar::xtile(ypc_w03, n = 5, wt = ponderador02_w03),
    quintil_w03 = factor(quintil_w03,
                         levels = 1:5,
                         labels = glue::glue('Quintil {1:5}')),
    quintil_w04 = statar::xtile(ypc_w04, n = 5, wt = ponderador02_w04),
    quintil_w04 = factor(quintil_w04,
                         levels = 1:5,
                         labels = glue::glue('Quintil {1:5}')),
    
    # Nivel educacional: 
    educ_w01 = factor(car::recode(m01_w01, recodes = "1:3 = 1; 4:5 = 2; 6:7 = 3; 8:10 = 4"),
                      levels = 1:4,
                      labels = c("Basica", "Media", "Tecnica", "Universitaria")),
    
    # Situacion ocupacional:
    ocup_w01 = factor(car::recode(m02_w01, "1:3 = 1; 6 = 2; c(4, 5, 8, 9) = 3; 7 = 4"),
                      levels = 1:4,
                      labels = c('Ocupado', 'Desempleado', 'Inactivo', 
                                 'Trabajo no remunerado')),
    ocup_w02 = factor(car::recode(m02_w02, "1:3 = 1; 6 = 2; c(4, 5, 8, 9) = 3; 7 = 4"),
                      levels = 1:4,
                      labels = c('Ocupado', 'Desempleado', 'Inactivo', 
                                 'Trabajo no remunerado')),
    ocup_w03 = factor(car::recode(m02_w03, "1:3 = 1; 6 = 2; c(4, 5, 8, 9) = 3; 7 = 4"),
                      levels = 1:4,
                      labels = c('Ocupado', 'Desempleado', 'Inactivo', 
                                 'Trabajo no remunerado')),
    ocup_w04 = factor(car::recode(m02_w04, "1:3 = 1; 6 = 2; c(4, 5, 8, 9) = 3; 7 = 4"),
                      levels = 1:4,
                      labels = c('Ocupado', 'Desempleado', 'Inactivo', 
                                 'Trabajo no remunerado')),
    
    
    # Sobrecarga de deudas
    sobrecarga_deuda_w01 = factor(car::recode(m43_w01, "1:2 = 1; 4:5 = 2; else = NA"),
                                  labels = c('Nada/No muy sobrecargado', 
                                             'Bastante/Muy sobrecargado')),
    
    sobrecarga_deuda_w03 = factor(car::recode(m43_w03, "1:2 = 1; 4:5 = 2; else = NA"),
                                  labels = c('Nada/No muy sobrecargado', 
                                             'Bastante/Muy sobrecargado')),
    
    # Apoyo social
    apoyo_social_w01 = factor(car::recode(s12_w01, "1:2 = 1; 3:4 = 2; 5 = 3"),
                              levels = 1:3,
                              labels = c('Apoyo social bajo', 'Apoyo social medio', 
                                         'Apoyo social alto')),
    apoyo_social_w03 = factor(car::recode(s12_w03, "1:2 = 1; 3:4 = 2; 5 = 3"),
                              levels = 1:3,
                              labels = c('Apoyo social bajo', 'Apoyo social medio', 
                                         'Apoyo social alto')),
    
    # Actividad fisica regular
    activ_fisica_w01 = factor(s04_w01 %in% 1:3,
                              labels = c('Sin actividad fisica regular', 
                                         'Con actividad fisica regular')),
    activ_fisica_w03 = factor(s04_w03 %in% 1:3,
                              labels = c('Sin actividad fisica regular', 
                                         'Con actividad fisica regular')),
    
    # Estado de salud
    salud_w01 = factor(car::recode(s03_w01, "1 = 1; 2:3 = 2; 4:5 = 3; else = NA"),
                       labels = c('Mala', 'Regular/buena', 'Muy buena/excelente')),
    salud_w02 = factor(car::recode(s03_w02, "1 = 1; 2:3 = 2; 4:5 = 3; else = NA"),
                       labels = c('Mala', 'Regular/buena', 'Muy buena/excelente')),
    salud_w03 = factor(car::recode(s03_w03, "1 = 1; 2:3 = 2; 4:5 = 3; else = NA"),
                       labels = c('Mala', 'Regular/buena', 'Muy buena/excelente')),
    salud_w04 = factor(car::recode(s03_w04, "1 = 1; 2:3 = 2; 4:5 = 3; else = NA"),
                       labels = c('Mala', 'Regular/buena', 'Muy buena/excelente')),
    salud_w05 = factor(car::recode(s03_w05, "1 = 1; 2:3 = 2; 4:5 = 3; else = NA"),
                       labels = c('Mala', 'Regular/buena', 'Muy buena/excelente')),
    
    # Fumador
    fuma_w01 = factor(s08_w01 >= 1,
                      labels = c('No fuma', 'Fumador')),
    fuma_w03 = factor(s08_w03 >= 1,
                      labels = c('No fuma', 'Fumador')),
    
    # IMC
    imc_w01 = ifelse(is_nsnr(s06_w01), NA, s06_w01) / (ifelse(is_nsnr(s05_w01), NA, s05_w01/100))^2,
    imc_w03 = ifelse(is_nsnr(s06_w03), NA, s06_w03) / (ifelse(is_nsnr(s05_w03), NA, s05_w03/100))^2,
    
    # Variables de personalidad
    
    s22_01_w02 = ifelse(is_nsnr(s22_01_w02), NA, s22_01_w02),
    s22_02_w02 = ifelse(is_nsnr(s22_02_w02), NA, s22_02_w02),
    s22_03_w02 = ifelse(is_nsnr(s22_03_w02), NA, s22_03_w02),
    s22_04_w02 = ifelse(is_nsnr(s22_04_w02), NA, s22_04_w02),
    s22_05_w02 = ifelse(is_nsnr(s22_05_w02), NA, s22_05_w02),
    s22_06_w02 = ifelse(is_nsnr(s22_06_w02), NA, s22_06_w02),
    s22_07_w02 = ifelse(is_nsnr(s22_07_w02), NA, s22_07_w02),
    s22_08_w02 = ifelse(is_nsnr(s22_08_w02), NA, s22_08_w02),
    s22_09_w02 = ifelse(is_nsnr(s22_09_w02), NA, s22_09_w02),
    s22_10_w02 = ifelse(is_nsnr(s22_10_w02), NA, s22_10_w02),
    
    extraversion_w02 = 6 - (s22_01_w02) + s22_06_w02 ,
    agreeableness_w02 = s22_02_w02 + 6 - s22_07_w02,
    conscientiousness_w02 = 6 - s22_03_w02 + s22_08_w02 ,
    neuroticism_w02 = 6 - s22_04_w02 + s22_09_w02,
    openness_w02 = 6 - s22_05_w02 + s22_10_w02,
    
    # Direccion vida
    proyecto_vida_w03 = factor(s30_01_w03 %in% 4:5,
                               labels = c('No decidido direccion de vida', 'Decidido direccion de vida')),
    
    # Desigualdad percibida
    
    percepcion_desigualdad_w01 = factor(c18_11_w01 %in% 4:5,
                                        labels = c('Desigualdad no muy grande', 
                                                   'Desigualdad muy grande')),
    percepcion_desigualdad_w02 = factor(c18_11_w02 %in% 4:5,
                                        labels = c('Desigualdad no muy grande', 
                                                   'Desigualdad muy grande')),
    percepcion_desigualdad_w03 = factor(c18_11_w03 %in% 4:5,
                                        labels = c('Desigualdad no muy grande', 
                                                   'Desigualdad muy grande')),
    percepcion_desigualdad_w04 = factor(c18_11_w04 %in% 4:5,
                                        labels = c('Desigualdad no muy grande', 
                                                   'Desigualdad muy grande')),
    
    # Estresores:
    
    estresor_01_w01 = s13_01_w01 == 1,
    estresor_02_w01 = s13_02_w01 == 1,
    estresor_03_w01 = s13_03_w01 == 1,
    estresor_04_w01 = s13_04_w01 == 1,
    estresor_06_w01 = s13_06_w01 == 1,
    
    estresor_01_w03 = s13_01_w03 == 1,
    estresor_02_w03 = s13_02_w03 == 1,
    estresor_03_w03 = s13_03_w03 == 1,
    estresor_04_w03 = s13_04_w03 == 1,
    estresor_06_w03 = s13_06_w03 == 1,
    
    estresor_n_w01 = estresor_01_w01 + estresor_02_w01 + estresor_03_w01 + estresor_04_w01 + estresor_06_w01,
    estresor_n_w03 = estresor_01_w03 + estresor_02_w03 + estresor_03_w03 + estresor_04_w03 + estresor_06_w03,
    
    estresor_n = estresor_n_w01 + estresor_n_w03
    
  )



####################################
### Datos ELSOC en formato LONG: ###
####################################

elsoc::load_elsoc(data = 'long')

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

save(elsoc_salud_wide, elsoc_salud_long, 
     file = file.path('..', 'inputs', 'datos_elsoc_preparados.RData'))


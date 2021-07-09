##############################################################################
# Codigo para la funcion de PAM
# 
# Gaston Quero - Sebastian Simondi                                
# 30/09/2020
################################################################################

getwd()

setwd("C:/Users/Usuario/OneDrive/Documentos/Tesis_Seba_Piazza")


# Paquetes 
library (lme4)
library (emmeans)
library ("car")
library ("nlmrt")
library ("easynls")
library ("ggplot2")       
library ("lattice")
library ("latticeExtra")
library (multcompView)
library (multcomp)
library ("dplyr")
library (ggjoy)
library ("ggridges")
library (hrbrthemes)
library(tidyverse)
library(forcats)
library("viridis")
library("lmerTest")
library(lubridate)
library(nycflights13)
library(nlme)
library(xtable)
library(stringr)
library(data.table)
library(svMisc)
library(ggpubr)
library("ggsci")
library("FactoMineR")
library("factoextra")
library("corrplot")


#C:\Users\Usuario\OneDrive\Documentos\Tesis_Seba_Piazza\Data\rawdata\preliminares\Resultados\Simondi 3\HHK1C3\curva
# Se carga los datos 

files.param <- dir(str_c("./Data/rawdata/preliminares/Resultados/parametros"), pattern = "*.asc")

pre.PAM.param <- bind_rows (lapply(files.param, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/preliminares/Resultados/parametros/", filt.raw) , 
                    delim =",", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = TRUE, 
                    col_types = list(col_character (), col_character(),col_character ()),
                    locale = default_locale(), 
                    na = "NA")
  dt <- dt %>%
       dplyr::mutate (pot=filt.raw)
  
}))

#pre.PAM.param$Parameter <- as.factor (pre.PAM.param$Parameter)

levels(pre.PAM.param$Parameter)

pre.PAM.param <- pre.PAM.param %>%
                 dplyr::mutate (Parameter = fct_recode (Parameter, "phi.psII.h"= "èPS2"))

param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.1 <- pre.PAM.param %>%
                   dplyr::filter (Parameter %in% param)%>%
                   dplyr::select (-X3)

head (pre.PAM.param.1)

write_delim(pre.PAM.param.1, path ="./Data/procdata/parametros_preliminares.txt" , 
            delim = ",", na = "NA")

### aca reingreso los datos 
PAM.param <- read_delim (file ="./Data/procdata/parametros_preliminares.txt", 
                  delim =",", quote = "\"",
                  escape_backslash = FALSE,
                  escape_double = TRUE, 
                  col_names = TRUE, 
                  col_types = NULL, 
                  locale = default_locale(), 
                  na = "NA")


PAM.param.1 <- PAM.param %>%
               dplyr::mutate (x = pot)%>%
               tidyr::separate(x , c("protocol", "genotype", "rep", NA), sep ="_")

unique (PAM.param.1$protocol)

PAM.param.1.s2 <- PAM.param.1 %>%
                  dplyr::filter (protocol ==  "Simondi2")


PAM.param.1.s3 <- PAM.param.1 %>%
                  dplyr::filter (protocol !=  "Simondi2")


pot.prueba <- PAM.param.1.s3 %>%
              dplyr::filter (pot == "Simondi3.1_HHK1C3_1_Param.asc") 

row.names(pot.prueba)

#### aca tiene que empezar la funcion 

pot.prueba <- pot.prueba %>%
              dplyr::mutate (ord = as.numeric(row.names (pot.prueba))) %>%
              dplyr::select (ord, everything())

head (pot.prueba)

t0 <- c(1,2)

Fo <- pot.prueba %>%
      dplyr::filter (ord %in% t0) %>%
      dplyr::summarise (Fo = mean (Value))

Fm <- pot.prueba %>%
      dplyr::filter (ord == 3) 

Fm.p <- pot.prueba %>%
           dplyr::filter (Parameter == "Fm'" )

Ft <- pot.prueba %>%
      dplyr::filter (Parameter == "Fs" ) 

quenching <- data.frame (Fm.prim = Fm.p$Value, Ft = Ft$Value)
  
phi.ps2 <- quenching  %>%
           dplyr::mutate (phiPS2 = (Fm.prim - Ft)/Fm.prim)

phi.NPQ <- quenching  %>%
           dplyr::mutate (phiNPQ = (Ft/Fm.prim) - Ft/Fm$Value)

phi.NO <- quenching  %>%
          dplyr::mutate (phiNO = Ft/Fm$Value)


quenching.analisis <- data.frame (phiPS2 = phi.ps2$phiPS2, 
                                  phiNPQ = phi.NPQ$phiNPQ,
                                  phiNO  = phi.NO$phiNO)

quenching.analisis <- quenching.analisis %>%
                      dplyr::mutate (sum.phi = phiPS2 + phiNPQ +  phiNO)

Fm.p.1 <- Fm.p %>%
          dplyr::filter ( ord == max(ord))

Ft.1 <- Ft %>%
        dplyr::filter ( ord == max(ord))


Fm.prim.prim <- pot.prueba %>%
                dplyr::filter ( ord > Fm.p.1$ord) %>%
                dplyr::filter (Parameter == "Fm")

relax <- data.frame ( Fm.pp = Fm.prim.prim$Value, 
                      Ft = Ft.1$Value, 
                      Fm.prim = Fm.p.1$Value)
800/60
phiNPQ.fast <- relax %>%
               dplyr::mutate (phiNPQfast = ((1/Fm.prim) - (1/Fm.pp))/(1/Ft))
   

phiNPQ.slow <- relax %>%
               dplyr::mutate (phiNPQslow = ((1/Fm.pp) - (1/Fm$Value))/(1/Ft))

relax.analisis <- data.frame (phiNPQfast = phiNPQ.fast$phiNPQfast,
                              phiNPQslow = phiNPQ.slow$phiNPQslow)
  
  
quantum.yield <- bind_cols (quenching.analisis,relax.analisis )







run.plot.curve <- function ( dt=NULL, sampling = NULL){
  
  dir.create (file.path ("Figures", "Plots.PAM"), showWarnings = FALSE)
  
  dt.1 <- dt %>%
    dplyr::select (-X3) %>%
    dplyr::rename (time= "Time (s)")%>%
    dplyr::mutate (time = dseconds (time))%>%
    dplyr::rename (fluorescence= "Data")
  
  PAM.curve <- ggscatter (dt.1 , x = "time", y = "fluorescence", 
                          title = str_c(sampling,unique(dt.1$pot),sep="_"), 
                          xlab = "time (s)",
                          ylab = "Fluorescence (a.u)",
                          point=FALSE) +
    geom_line(color = "black", size = 0.5)
  
  ggexport (PAM.curve, filename = str_c("./Figures/Plots.PAM/",sampling,unique(dt.1$pot), ".tiff"),
            width = 700, height = 500)

  print (PAM.curve)          
}

X <- lapply(files.HHK1C3.S3, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves.HHK1C3.s3 %>%
         dplyr::filter (pot == filt.pot)
         
  run.plot.curve (dt=pot,sampling = "HHK1C3.s3" )    
  
  
})



files.GT.1.S3.1 <- dir(str_c("./Data/rawdata/GT.1/simondi_3.1/curva_s3.1"), pattern = "*.asc")


PAM.curves.GT.1.s3.1 <- bind_rows (lapply(files.GT.1.S3.1, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/GT.1/simondi_3.1/curva_s3.1/", filt.raw) , 
                    delim =",", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = TRUE, 
                    col_types = NULL,
                    locale = default_locale(), 
                    na = "NA")
  dt <- dt %>%
    dplyr::mutate (pot=filt.raw)
  
}))


X.1 <- lapply(files.GT.1.S3.1, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves.GT.1.s3.1 %>%
    dplyr::filter (pot == filt.pot)
  
  run.plot.curve (dt=pot,sampling = "GT.1.s3.1" )    
  
  
})




















files.diciembre <- dir(str_c("./Data/rawdata/PAM.diciembre/curva"), pattern = "*.asc")
PAM.curves.diciembre <- bind_rows (lapply(files.diciembre, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/PAM.diciembre/curva/", filt.raw) , 
                    delim =",", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = TRUE, 
                    col_types = NULL,
                    locale = default_locale(), 
                    na = "NA")
  dt <- dt %>%
    dplyr::mutate (pot=filt.raw)
  
}))


Y <- lapply(files.diciembre, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves.diciembre %>%
    dplyr::filter (pot == filt.pot)
  
  run.plot.curve (dt=pot,sampling = "diciembre" )    
  
  
})


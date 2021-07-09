##############################################################################
# Codigo para la figura de protocolos de PAM 
# 
# Gaston Quero - Sebastian Simondi                                
# 25/03/2021
################################################################################

getwd()

setwd ("C:/Users/Usuario/OneDrive/Documentos/Tesis_Seba_Piazza")

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

### Aca reingreso los de las curvas de los protocolos ####

dir.data.base <- ("C:/Users/Usuario/OneDrive/Documentos/Paper_fotosintesis_dinámica")

####### piazza 400 #####

PAM.curves.400.2 <- read_delim (file = str_c(dir.data.base, "/Data/procdata/SP_1/PAM.curves.SP.1.pz400.2.txt") , 
                                delim =",", quote = "\"",
                                escape_backslash = FALSE,
                                escape_double = TRUE, 
                                col_names = TRUE, 
                                col_types = NULL,
                                locale = default_locale(), 
                                na = "NA")

PAM.curves.400 <- read_delim (file = str_c(dir.data.base, "/Data/procdata/SP_1/PAM.curves.SP.1.pz400.txt") , 
                                delim =",", quote = "\"",
                                escape_backslash = FALSE,
                                escape_double = TRUE, 
                                col_names = TRUE, 
                                col_types = NULL,
                                locale = default_locale(), 
                                na = "NA")


PAM.curves.800.2 <- read_delim (file = str_c(dir.data.base, "/Data/procdata/SP_1/PAM.curves.SP.1.pz800.2.txt") , 
                              delim =",", quote = "\"",
                              escape_backslash = FALSE,
                              escape_double = TRUE, 
                              col_names = TRUE, 
                              col_types = NULL,
                              locale = default_locale(), 
                              na = "NA")


PAM.curves.800 <- read_delim (file = str_c(dir.data.base, "/Data/procdata/SP_1/PAM.curves.SP.1.pz800.txt") , 
                                delim =",", quote = "\"",
                                escape_backslash = FALSE,
                                escape_double = TRUE, 
                                col_names = TRUE, 
                                col_types = NULL,
                                locale = default_locale(), 
                                na = "NA")


### Inv ###

PAM.curves.400.2.inv <- PAM.curves.400.2 %>%
                        dplyr::filter (ambiente == "Inv")


PAM.curves.400.inv <- PAM.curves.400 %>%
                      dplyr::filter (ambiente == "Inv")


PAM.curves.800.2.inv <- PAM.curves.800.2 %>%
                        dplyr::filter (ambiente == "Inv")


PAM.curves.800.inv <- PAM.curves.800 %>%
                      dplyr::filter (ambiente == "Inv")


###
Inv <- bind_rows (PAM.curves.400.2.inv,
                         PAM.curves.400.inv,
                         PAM.curves.800.2.inv,
                         PAM.curves.800.inv)
 
## plot de la firma del PAM
dt.1  <- Inv %>%
                 dplyr::filter (genotype == "HHK1C3") %>%
                 dplyr::filter (rep == 4)



dt.1 <- dt.1 %>%
        dplyr::rename (time= "Time (s)")%>%
        dplyr::mutate (time = dseconds (time))%>%
        dplyr::rename (fluorescence= "Data")


dt.fm <- dt.1 %>%
         dplyr::filter ( fluorescence == max (fluorescence))

t.fvfm <- dt.fm$time [1] + 0

dt.act <- dt.1 %>%
         dplyr::filter (  time >= t.fvfm + 3.2) %>%
         dplyr::filter ( fluorescence == max (fluorescence))

t.act <- 43.23 + 0


t.c1 <- t.act + 60 +1.6
t.c2 <- t.c1 + 15 + 2
t.c3 <- t.c2 + 15 + 3
t.c4 <- t.c3 + 15 + 3
t.c5 <- t.c4 + 15 + 3
t.c6 <- t.c5 + 15 + 2
t.c7 <- t.c6 + 15 + 2
t.c8 <- t.c7 + 15 + 2
t.c9 <- t.c8 + 15 + 3
t.c10 <- t.c9 + 15  + 2
t.c11 <- t.c10 + 15 + 3
t.c12 <- t.c11 + 15 + 2
t.c13 <- t.c12 + 15 + 3
t.c14 <- t.c13 + 15 + 2
t.c15 <- t.c14 + 15 + 2
t.c16 <- t.c15 + 15 + 2
t.c17 <- t.c16 + 15 + 2
t.c18 <- t.c17 + 15 + 2
t.c19 <- t.c18 + 15 + 2
t.c20 <- t.c19 + 15 
t.c21 <- t.c20 + 15 + 4
t.c22 <- t.c21 + 30 + 4
t.c23 <- t.c22 + 30 + 4
t.c24 <- t.c23 + 30 + 4
t.c25 <- t.c24 + 30 + 5
t.c26 <- t.c25 + 30 + 5
t.c27 <- t.c26 + 30 + 5
t.c28 <- t.c27 + 30 + 5
t.c29 <- t.c28 + 30 + 5
t.c30 <- t.c29 + 30 + 5


PAM.curve <- ggscatter (dt.1 , x = "time", y = "fluorescence", facet.by = "ppfd.act",
                         #title = str_c(sampling,protocolo ,unique(dt.1$pot),sep="_"), 
                         xlab = "time (s)",
                         ylab = "Fluorescence (a.u)",
                         point=FALSE) +
           geom_line  (color = "black", size = 0.5)




ggexport (PAM.curve, filename = str_c("./Figures/","_", unique(dt.1$pot), ".tiff"),
          width = 700, height = 500)



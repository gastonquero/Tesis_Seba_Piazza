##########################################################
# Codigo para el analisis de datos de Sebastian Piazza   #
#                                                        #
# Gaston Quero - Sebastian Simondi                       #          
# 1/12/2020                                             #
##########################################################

###########################
# Los datos fueron enviados 
# de:	Sebastián Piazza <spiazzajorcin@gmail.com>
#  para:	Gastón Quero <gastonquero@gmail.com>
#  fecha:	25 nov 2020 12:40
# asunto:	Resultados medición
# enviado por:	gmail.com


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

###
#C:\Users\Usuario\OneDrive\Documentos\Tesis_Seba_Piazza\Data\rawdata\preliminares\Resultados\Simondi 3\HHK1C3\curva
# Se carga los datos 

files.curva.400 <- dir(str_c("./Data/rawdata/Piazza_400/curvas_pza400"), pattern = "*.asc")

PAM.curves.400 <- bind_rows (lapply(files.curva.400, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/Piazza_400/curvas_pza400/", filt.raw) , 
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
  

lapply(files.curva.400, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves.400 %>%
         dplyr::filter (pot == filt.pot)
         
  run.plot.curve (dt=pot, sampling = "piazza.400" )    
  
  
})

unique (PAM.curves.400$pot)
head (PAM.curves.400)

PAM.curves.400  <- PAM.curves.400 %>%
                   dplyr::mutate (x = pot)%>%
                   tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") %>%
                   dplyr::select (-X3)

PAM.curves.400$ambiente


Inv.400 <- PAM.curves.400 %>%
           dplyr::filter ( ambiente == "Inv") %>%
           dplyr::mutate ( ppfd = "400" ) %>%
           dplyr::select (ambiente,  ppfd, genotype, rep, pot,`Time (s)`,  Data )

head (LB.400)
LB <-  PAM.curves.400 %>%
       dplyr::filter (ambiente != "Inv")

LB1 <- LB %>% 
      separate(ppfd, c("ppfd", NA))%>%
      dplyr::select (ambiente,  ppfd, genotype, rep, pot,`Time (s)`,  Data )

PAM.curves.400.1 <- bind_rows (Inv.400, LB1 )
  
PAM.curves.400.1 <- PAM.curves.400.1 %>%
                    dplyr::mutate (ppfd.act = 4.03e2)
 
                
write_delim (PAM.curves.400.1, path ="./Data/procdata/PAM.curves.piazza.400.txt" , 
             delim = ",", na = "NA")



## Aca cargo los datos de los parametros
files.param.400 <- dir(str_c("./Data/rawdata/Piazza_400/parametros_pza400"), pattern = "*.asc")

pre.PAM.param.400 <- bind_rows (lapply(files.param.400, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/Piazza_400/parametros_pza400/", filt.raw) , 
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



param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.400.1 <- pre.PAM.param.400 %>%
                       dplyr::filter (Parameter %in% param)%>%
                       dplyr::select (-X3)

head (pre.PAM.param.400.1)

write_delim (pre.PAM.param.400.1, file ="./Data/procdata/pre.PAM.param.400.1.txt" , 
            delim = ",", na = "NA")

### aca reingreso los datos 
PAM.param.400 <- read_delim (file ="./Data/procdata/pre.PAM.param.400.1.txt", 
                         delim =",", quote = "\"",
                         escape_backslash = FALSE,
                         escape_double = TRUE, 
                         col_names = TRUE, 
                         col_types = NULL, 
                         locale = default_locale(), 
                         na = "NA")

unique ( PAM.param.400$pot)

PAM.param.400.1 <- PAM.param.400  %>%
                   dplyr::mutate (x = pot)%>%
                   tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") 


par.Inv.400 <- PAM.param.400.1 %>%
               dplyr::filter ( ambiente == "Inv") %>%
               dplyr::mutate ( ppfd = "400" ) %>%
               dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter, Value )


par.LB.400 <-  PAM.param.400.1  %>%
               dplyr::filter (ambiente != "Inv")

par.LB1.400 <- par.LB.400 %>% 
               separate(ppfd, c("ppfd", NA))%>%
               dplyr::select (ambiente,  ppfd, genotype, rep, pot,Parameter, Value)

PAM.param.400.1a <- bind_rows (par.Inv.400, par.LB1.400 )

PAM.param.400.1a <- PAM.param.400.1a %>%
                    dplyr::mutate (ppfd.act = 4.03e2)


write_delim (PAM.param.400.1a, file ="./Data/procdata/PAM.param.400.1a.txt" , 
             delim = ",", na = "NA")


################ piazza 400.2 
files.curva.400.2 <- dir(str_c("./Data/rawdata/Piazza_400.2/curvas_pza400.2"), pattern = "*.asc")

PAM.curves.400.2 <- bind_rows (lapply(files.curva.400.2, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/Piazza_400.2/curvas_pza400.2/", filt.raw) , 
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


lapply(files.curva.400.2, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves.400.2 %>%
    dplyr::filter (pot == filt.pot)
  
  run.plot.curve (dt=pot, sampling = "piazza.400.2" )    
  
  
})

unique (PAM.curves.400.2$pot)
head (PAM.curves.400.2)

PAM.curves.400.2  <- PAM.curves.400.2 %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") %>%
  dplyr::select (-X3)

PAM.curves.400.2$ambiente


Inv.400.2 <- PAM.curves.400.2 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "400.2" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,`Time (s)`,  Data )

head (LB.400.2)
LB <-  PAM.curves.400.2 %>%
  dplyr::filter (ambiente != "Inv")

LB1 <- LB %>% 
  separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,`Time (s)`,  Data )

PAM.curves.400.2.1 <- bind_rows (Inv.400.2, LB1 )

PAM.curves.400.2.1 <- PAM.curves.400.2.1 %>%
                      dplyr::mutate (ppfd.act = 200)


write_delim (PAM.curves.400.2.1, file ="./Data/procdata/PAM.curves.piazza.400.2.txt" , 
             delim = ",", na = "NA")



## Aca cargo los datos de los parametros
files.param.400.2 <- dir(str_c("./Data/rawdata/Piazza_400.2/parametros_pza400.2"), pattern = "*.asc")

pre.PAM.param.400.2 <- bind_rows (lapply(files.param.400.2, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/Piazza_400.2/parametros_pza400.2/", filt.raw) , 
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



param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.400.2.1 <- pre.PAM.param.400.2 %>%
  dplyr::filter (Parameter %in% param)%>%
  dplyr::select (-X3)

head (pre.PAM.param.400.2.1)

write_delim (pre.PAM.param.400.2.1, file ="./Data/procdata/pre.PAM.param.400.2.1.txt" , 
             delim = ",", na = "NA")

### aca reingreso los datos 
PAM.param.400.2 <- read_delim (file ="./Data/procdata/pre.PAM.param.400.2.1.txt", 
                             delim =",", quote = "\"",
                             escape_backslash = FALSE,
                             escape_double = TRUE, 
                             col_names = TRUE, 
                             col_types = NULL, 
                             locale = default_locale(), 
                             na = "NA")

unique ( PAM.param.400.2$pot)

PAM.param.400.2.1 <- PAM.param.400.2  %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") 


par.Inv.400.2 <- PAM.param.400.2.1 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "400.2" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter, Value )


par.LB.400.2 <-  PAM.param.400.2.1  %>%
  dplyr::filter (ambiente != "Inv")

par.LB1.400.2 <- par.LB.400.2 %>% 
  separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,Parameter, Value)

PAM.param.400.2.1a <- bind_rows (par.Inv.400.2, par.LB1.400.2 )

PAM.param.400.2.1a <- PAM.param.400.2.1a %>%
                      dplyr::mutate (ppfd.act = 200)


write_delim (PAM.param.400.2.1a, file ="./Data/procdata/PAM.param.400.2.1a.txt" , 
             delim = ",", na = "NA")

############## 800 ### 
files.curva.800 <- dir(str_c("./Data/rawdata/Piazza_800/curvas_pza800"), pattern = "*.asc")

PAM.curves.800 <- bind_rows (lapply(files.curva.800, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/Piazza_800/curvas_pza800/", filt.raw) , 
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


lapply(files.curva.800, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves.800 %>%
    dplyr::filter (pot == filt.pot)
  
  run.plot.curve (dt=pot, sampling = "piazza.800" )    
  
  
})

unique (PAM.curves.800$pot)
head (PAM.curves.800)

PAM.curves.800  <- PAM.curves.800 %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") %>%
  dplyr::select (-X3)

PAM.curves.800$ambiente


Inv.800 <- PAM.curves.800 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "400" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,`Time (s)`,  Data )


LB <-  PAM.curves.800 %>%
  dplyr::filter (ambiente != "Inv")

LB1 <- LB %>% 
  separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,`Time (s)`,  Data )

PAM.curves.800.1 <- bind_rows (Inv.800, LB1 )

PAM.curves.800.1 <- PAM.curves.800.1 %>%
  dplyr::mutate (ppfd.act = 1.7e3)


write_delim (PAM.curves.800.1, path ="./Data/procdata/PAM.curves.piazza.800.txt" , 
             delim = ",", na = "NA")



## Aca cargo los datos de los parametros
files.param.800 <- dir(str_c("./Data/rawdata/Piazza_800/parametros_pza800"), pattern = "*.asc")

pre.PAM.param.800 <- bind_rows (lapply(files.param.800, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/Piazza_800/parametros_pza800/", filt.raw) , 
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



param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.800.1 <- pre.PAM.param.800 %>%
  dplyr::filter (Parameter %in% param)%>%
  dplyr::select (-X3)

head (pre.PAM.param.800.1)

write_delim (pre.PAM.param.800.1, file ="./Data/procdata/pre.PAM.param.800.1.txt" , 
             delim = ",", na = "NA")

### aca reingreso los datos 
PAM.param.800 <- read_delim (file ="./Data/procdata/pre.PAM.param.800.1.txt", 
                             delim =",", quote = "\"",
                             escape_backslash = FALSE,
                             escape_double = TRUE, 
                             col_names = TRUE, 
                             col_types = NULL, 
                             locale = default_locale(), 
                             na = "NA")

unique ( PAM.param.800$pot)

PAM.param.800.1 <- PAM.param.800  %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") 


par.Inv.800 <- PAM.param.800.1 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "400" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter, Value )


par.LB.800 <-  PAM.param.800.1  %>%
  dplyr::filter (ambiente != "Inv")

par.LB1.800 <- par.LB.800 %>% 
  separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,Parameter, Value)

PAM.param.800.1a <- bind_rows (par.Inv.800, par.LB1.800 )

PAM.param.800.1a <- PAM.param.800.1a %>%
  dplyr::mutate (ppfd.act = 1.7e3)


write_delim (PAM.param.800.1a, file ="./Data/procdata/PAM.param.800.1a.txt" , 
             delim = ",", na = "NA")


################ piazza 800.2 
files.curva.800.2 <- dir(str_c("./Data/rawdata/Piazza_800.2/curvas_pza800.2"), pattern = "*.asc")

PAM.curves.800.2 <- bind_rows (lapply(files.curva.800.2, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/Piazza_800.2/curvas_pza800.2/", filt.raw) , 
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


lapply(files.curva.800.2, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves.800.2 %>%
    dplyr::filter (pot == filt.pot)
  
  run.plot.curve (dt=pot, sampling = "piazza.800.2" )    
  
  
})

unique (PAM.curves.800.2$pot)
head (PAM.curves.800.2)

PAM.curves.800.2  <- PAM.curves.800.2 %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") %>%
  dplyr::select (-X3)

PAM.curves.800.2$ambiente


Inv.800.2 <- PAM.curves.800.2 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "400" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,`Time (s)`,  Data )

head (LB.800.2)
LB <-  PAM.curves.800.2 %>%
  dplyr::filter (ambiente != "Inv")

LB1 <- LB %>% 
  separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,`Time (s)`,  Data )

PAM.curves.800.2.1 <- bind_rows (Inv.800.2, LB1 )

PAM.curves.800.2.1 <- PAM.curves.800.2.1 %>%
  dplyr::mutate (ppfd.act = 850)


write_delim (PAM.curves.800.2.1, file ="./Data/procdata/PAM.curves.piazza.800.2.txt" , 
             delim = ",", na = "NA")



## Aca cargo los datos de los parametros
files.param.800.2 <- dir(str_c("./Data/rawdata/Piazza_800.2/parametros_pza800.2"), pattern = "*.asc")

pre.PAM.param.800.2 <- bind_rows (lapply(files.param.800.2, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/Piazza_800.2/parametros_pza800.2/", filt.raw) , 
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



param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.800.2.1 <- pre.PAM.param.800.2 %>%
  dplyr::filter (Parameter %in% param)%>%
  dplyr::select (-X3)

head (pre.PAM.param.800.2.1)

write_delim (pre.PAM.param.800.2.1, file ="./Data/procdata/pre.PAM.param.800.2.1.txt" , 
             delim = ",", na = "NA")

### aca reingreso los datos 
PAM.param.800.2 <- read_delim (file ="./Data/procdata/pre.PAM.param.800.2.1.txt", 
                               delim =",", quote = "\"",
                               escape_backslash = FALSE,
                               escape_double = TRUE, 
                               col_names = TRUE, 
                               col_types = NULL, 
                               locale = default_locale(), 
                               na = "NA")

unique ( PAM.param.800.2$pot)

PAM.param.800.2.1 <- PAM.param.800.2  %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") 


par.Inv.800.2 <- PAM.param.800.2.1 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "800.2" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter, Value )


par.LB.800.2 <-  PAM.param.800.2.1  %>%
  dplyr::filter (ambiente != "Inv")

par.LB1.800.2 <- par.LB.800.2 %>% 
  separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot,Parameter, Value)

PAM.param.800.2.1a <- bind_rows (par.Inv.800.2, par.LB1.800.2 )

PAM.param.800.2.1a <- PAM.param.800.2.1a %>%
  dplyr::mutate (ppfd.act = 850)


write_delim (PAM.param.800.2.1a, file ="./Data/procdata/PAM.param.800.2.1a.txt" , 
             delim = ",", na = "NA")








######
### hay que ir por protocolo
sampling = "preliminar_"

#param.s3 <- PAM.param.1 %>%
            #dplyr::filter ( protocol == "Simondi3")

list.pot.s3 <- unique (param.s3$pot)
unique (PAM.param.1$protocol)

dt = PAM.param.1
protocol = "Simondi3.1"


run_quenching_analisis <-  function ( dt = NULL, sampling = NULL, protocol = NULL) {
  
  prt <- protocol
  param.protocol <- dt %>%
                    dplyr::filter ( protocol == prt)
  
  list.pot <- unique (param.protocol$pot)
  
  df.pot <- bind_rows (lapply( list.pot, function (filt.pot){
    
    #filt.pot = "Simondi3.1_HHK1C3_1_Param.asc"
    
    pot.par <- dt %>%
               dplyr::filter (pot == filt.pot )
    
    pot.par <- pot.par %>%
               dplyr::mutate (ord = as.numeric(row.names (pot.par))) %>%
               dplyr::select (ord, everything())
    
    id <- pot.par %>%
          dplyr::select (c(pot, protocol, genotype, rep )) %>%
          dplyr::select (protocol, genotype ,rep,  pot  )
    ## parametros inciales
    fluor.min <- pot.par %>%
      dplyr::filter (ord < 3) %>%
      dplyr::summarise (Fo = mean (Value))
    
    fluor.max <- pot.par %>%
      dplyr::filter (ord == 3) %>%
      dplyr::mutate (Fm = Value)%>%
      dplyr::select (Fm)
    
    fluor.max.prim <- pot.par %>%
      dplyr::filter (Parameter == "Fm'" ) %>%
      dplyr::filter (ord != max(ord))
    
    fluor <- pot.par  %>%
      dplyr::filter (Parameter == "Fs") %>%
      dplyr::filter (ord != max(ord))
    
    phi.ps2 <- data.frame ( phiPS2 = (fluor.max.prim$Value - fluor$Value)/fluor.max.prim$Value)
    
    phi.NPQ <- data.frame ( phiNPQ = (fluor$Value/fluor.max.prim$Value) -(fluor$Value/fluor.max$Fm))
    
    phi.NO <- data.frame ( phiNO =fluor$Value/fluor.max$Fm)
    
    quenching.analisis <- bind_cols(phi.ps2, phi.NPQ, phi.NO )
    
    quenching.analisis <- quenching.analisis %>%
                          dplyr::mutate  (ciclo = row.names (quenching.analisis)) 
    
    x <- length(unique (quenching.analisis$ciclo))-1
    x1 <- c(1:x * 30) + 60
    time <- c(60,x1 )
    
    quenching.analisis <- quenching.analisis %>%
                           dplyr::mutate (time.min = time/60)
    
    quenching.analisis.1 <- quenching.analisis %>%
                             pivot_longer(-c(ciclo, time.min), names_to = "quantum.yield" ,values_to = "phi")
    
    quenching.analisis.1 <- quenching.analisis.1  %>%
                            dplyr::mutate(quantum.yield = factor(quantum.yield, levels = c("phiPS2","phiNPQ","phiNO")))%>%
                            dplyr::mutate (phi = round (phi,2))
    
    #ggbarplot(quenching.analisis.1 , "ciclo", "phi",title = unique (id$pot),
    # fill = "quantum.yield", color = "quantum.yield", 
    #palette = c("gray18","gray38" ,"gray58"),
    #label = TRUE, lab.col = "white", lab.pos = "in")
    
    phi.plot <- ggbarplot(quenching.analisis.1 , "time.min", "phi",title =str_c(sampling,unique (id$pot), sep="_"),
                          ylab="time (min)",
                          fill = "quantum.yield", color = "quantum.yield", 
                          palette = c("gray18","gray38" ,"gray58"),
                          label = TRUE, lab.col = "white", lab.pos = "in")
    
    ggexport (phi.plot, filename = str_c("./Figures/Plots.PAM/",sampling,"phi.plot", unique(id$pot), ".tiff"),
              width = 700, height = 500)
    
    print (phi.plot)
    xid <-id [1,] 
    
    quenching.analisis.2 <- cbind (quenching.analisis, xid , sampling = sampling)
    
    return (quenching.analisis.2 )
  }))
  
  
  return (df.pot)
}


preliminar_S3 <- run_quenching_analisis ( dt = PAM.param.1, sampling ="preliminar",protocol = "Simondi3")

preliminar_S3.1 <- run_quenching_analisis ( dt = PAM.param.1, sampling ="preliminar",protocol = "Simondi3.1")

preliminar_S3.2 <- run_quenching_analisis ( dt = PAM.param.1, sampling ="preliminar",protocol = "Simondi3.2")


preliminar_relax_S3 <- run_relax_analisis ( dt = PAM.param.1, sampling = "preliminar", protocol = "Simondi3", rep =TRUE)



preliminar_relax_S3.1 <- run_relax_analisis ( dt = PAM.param.1, sampling = "preliminar", protocol = "Simondi3.1", rep =TRUE)


preliminar_relax_S3.2 <- run_relax_analisis ( dt = PAM.param.1, sampling = "preliminar", protocol = "Simondi3.2", rep =TRUE)












pot.par <- PAM.param.1 %>%
             dplyr::filter (pot == "Simondi3_HHK1C3_1_Param.asc")
  
  
#}
plot.phi.s3 

pot.par <- PAM.param.1 %>%
           dplyr::filter (pot == "Simondi3_HHK1C3_1_Param.asc")

pot.curv <- PAM.curves.preliminar %>%
            dplyr::filter (protocol == unique (pot.par$protocol)) %>%
            dplyr::filter (genotype == unique (pot.par$genotype)) %>%
            dplyr::filter (rep == unique (pot.par$rep))



pot.curv.1 <- pot.curv  %>%
              dplyr::rename (time= "Time (s)")%>%
              dplyr::mutate (time = dseconds (time))%>%
              dplyr::rename (fluorescence= "Data")

plot.curve <- ggscatter (pot.curv.1 , x = "time", y = "fluorescence", 
                        title = str_c(sampling,unique(pot.curv.1$pot),sep="_"), 
                        xlab = "time (s)",
                        ylab = "Fluorescence (a.u)",
                        point=FALSE) +
              geom_line(color = "black", size = 0.5) +
              geom_vline(xintercept = 100, color = "black", size = 0.5)


print (plot.curve)      

# s3   = actinica 30
# s3.1 = actinica 45



pot.par <- pot.par %>%
           dplyr::mutate (ord = as.numeric(row.names (pot.par))) %>%
           dplyr::select (ord, everything())

id <- pot.par %>%
      dplyr::select (c(pot, protocol, genotype, rep )) %>%
      dplyr::select (protocol, genotype ,rep,  pot  )


## parametros inciales
fluor.min <- pot.par %>%
             dplyr::filter (ord < 3) %>%
             dplyr::summarise (Fo = mean (Value))

fluor.max <- pot.par %>%
             dplyr::filter (ord == 3) %>%
             dplyr::mutate (Fm = Value)%>%
             dplyr::select (Fm)

fluor.max.prim <- pot.par %>%
                  dplyr::filter (Parameter == "Fm'" ) %>%
                  dplyr::filter (ord != max(ord))

fluor <- pot.par  %>%
         dplyr::filter (Parameter == "Fs") %>%
         dplyr::filter (ord != max(ord))

phi.ps2 <- data.frame ( phiPS2 = (fluor.max.prim$Value - fluor$Value)/fluor.max.prim$Value)

phi.NPQ <- data.frame ( phiNPQ = (fluor$Value/fluor.max.prim$Value) -(fluor$Value/fluor.max$Fm))

phi.NO <- data.frame ( phiNO =fluor$Value/fluor.max$Fm)

quenching.analisis <- bind_cols(phi.ps2, phi.NPQ, phi.NO )

x <- length(unique (quenching.analisis.1$ciclo))-1
x1 <- c(1:x * 30) + 60
time <- c(60,x1 )

quenching.analisis <- quenching.analisis %>%
                      dplyr::mutate  (ciclo = row.names (quenching.analisis)) %>%
                      dplyr::mutate (time = time/60)

quenching.analisis.1 <- quenching.analisis %>%
                        pivot_longer(-c(ciclo, time), names_to = "quantum.yield" ,values_to = "phi")

quenching.analisis.1 <- quenching.analisis.1  %>%
                        dplyr::mutate(quantum.yield = factor(quantum.yield, levels = c("phiPS2","phiNPQ","phiNO")))%>%
                        dplyr::mutate (phi = round (phi,2))

#ggbarplot(quenching.analisis.1 , "ciclo", "phi",title = unique (id$pot),
         # fill = "quantum.yield", color = "quantum.yield", 
          #palette = c("gray18","gray38" ,"gray58"),
          #label = TRUE, lab.col = "white", lab.pos = "in")

phi.plot <- ggbarplot(quenching.analisis.1 , "time", "phi",title = unique (id$pot),
          fill = "quantum.yield", color = "quantum.yield", 
          palette = c("gray18","gray38" ,"gray58"),
          label = TRUE, lab.col = "white", lab.pos = "in")

ggexport (phi.plot, filename = str_c("./Figures/Plots.PAM/",sampling,"phi.plot", unique(id$pot), ".tiff"),
          width = 700, height = 500)

print (phi.plot)
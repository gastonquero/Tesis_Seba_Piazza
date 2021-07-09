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
#################################3


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

### Piazza.400

files.param.pza400 <- dir(str_c("./Data/rawdata/Piazza_400/parametros_pza400"), pattern = "*.asc")

pre.PAM.param.pza400  <- bind_rows (lapply(files.param.pza400, function (filt.raw) {
  
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

#pre.PAM.param$Parameter <- as.factor (pre.PAM.param$Parameter)

head (pre.PAM.param.pza400  )

param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.pza400.1 <- pre.PAM.param.pza400  %>%
                           dplyr::filter (Parameter %in% param)%>%
                           dplyr::select (-X3)


write_delim (pre.PAM.param.pza400.1 , file ="./Data/procdata/pre.PAM.param.pza400.1.txt" , 
            delim = ",", na = "NA")



## Aca cargo los datos de los parametros
#C:\Users\Usuario\OneDrive\Documentos\Tesis_Seba_Piazza\Data\rawdata\Piazza_400\parametros_pza400

### Piazza.400

files.param.pza400.2 <- dir(str_c("./Data/rawdata/Piazza_400.2/parametros_pza400.2"), pattern = "*.asc")

pre.PAM.param.pza400.2  <- bind_rows (lapply(files.param.pza400.2, function (filt.raw) {
  
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

#pre.PAM.param$Parameter <- as.factor (pre.PAM.param$Parameter)

head (pre.PAM.param.pza400.2  )

param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.pza400.2.1 <- pre.PAM.param.pza400.2  %>%
  dplyr::filter (Parameter %in% param)%>%
  dplyr::select (-X3)


write_delim (pre.PAM.param.pza400.2.1 , file ="./Data/procdata/pre.PAM.param.pza400.2.1.txt" , 
             delim = ",", na = "NA")



### Piazza.800

files.param.pza800 <- dir(str_c("./Data/rawdata/Piazza_800/parametros_pza800"), pattern = "*.asc")

pre.PAM.param.pza800  <- bind_rows (lapply(files.param.pza800, function (filt.raw) {
  
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

#pre.PAM.param$Parameter <- as.factor (pre.PAM.param$Parameter)

head (pre.PAM.param.pza800  )

param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.pza800.1 <- pre.PAM.param.pza800  %>%
  dplyr::filter (Parameter %in% param)%>%
  dplyr::select (-X3)


write_delim (pre.PAM.param.pza800.1 , file ="./Data/procdata/pre.PAM.param.pza800.1.txt" , 
             delim = ",", na = "NA")


### Piazza.800.2

files.param.pza800.2 <- dir(str_c("./Data/rawdata/Piazza_800.2/parametros_pza800.2"), pattern = "*.asc")

pre.PAM.param.pza800.2  <- bind_rows (lapply(files.param.pza800.2, function (filt.raw) {
  
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

#pre.PAM.param$Parameter <- as.factor (pre.PAM.param$Parameter)

head (pre.PAM.param.pza800.2  )

param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.pza800.2.1 <- pre.PAM.param.pza800.2  %>%
  dplyr::filter (Parameter %in% param)%>%
  dplyr::select (-X3)


write_delim (pre.PAM.param.pza800.2.1 , file ="./Data/procdata/pre.PAM.param.pza800.2.1.txt" , 
             delim = ",", na = "NA")



### aca reingreso los datos 
PAM.param.pza400 <- read_delim (file ="./Data/procdata/pre.PAM.param.pza400.1.txt", 
                         delim =",", quote = "\"",
                         escape_backslash = FALSE,
                         escape_double = TRUE, 
                         col_names = TRUE, 
                         col_types = NULL, 
                         locale = default_locale(), 
                         na = "NA")


PAM.param.pza400.1 <- PAM.param.pza400 %>%
                    dplyr::mutate (x = pot)%>%
                    tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") 

head (PAM.param.pza400.1)

Inv.param.400 <- PAM.param.pza400.1 %>%
                 dplyr::filter ( ambiente == "Inv") %>%
                 dplyr::mutate ( ppfd = "400" ) %>%
                 dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter ,  Value )


LB.param.400  <-  PAM.param.pza400.1  %>%
                  dplyr::filter (ambiente != "Inv")

LB.param.400  <- LB.param.400  %>% 
                 tidyr::separate(ppfd, c("ppfd", NA))%>%
                  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter , Value )


PAM.piazza.400.1 <- bind_rows (Inv.param.400, LB.param.400 )

PAM.piazza.400.1  <- PAM.piazza.400.1 %>%
                     dplyr::mutate (ppfd.act = 4.03e2)


PAM.piazza.400.1 <- PAM.piazza.400.1 %>%
                    dplyr::rename (sampling = ambiente) %>%
                    dplyr::mutate (protocol = "piazza.400" )


write_delim (PAM.piazza.400.1, path ="./Data/procdata/PAM.piazza.400.1.txt" , 
              delim = ",", na = "NA")

#### Piazza 400.2

PAM.param.pza400.2 <- read_delim (file ="./Data/procdata/pre.PAM.param.pza400.2.1.txt", 
                                delim =",", quote = "\"",
                                escape_backslash = FALSE,
                                escape_double = TRUE, 
                                col_names = TRUE, 
                                col_types = NULL, 
                                locale = default_locale(), 
                                na = "NA")


PAM.param.pza400.2.1 <- PAM.param.pza400.2 %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") 

head (PAM.param.pza400.2.1)

Inv.param.400.2 <- PAM.param.pza400.2.1 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "400.2" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter ,  Value )


LB.param.400.2  <-  PAM.param.pza400.2.1  %>%
  dplyr::filter (ambiente != "Inv")

LB.param.400.2  <- LB.param.400.2  %>% 
  tidyr::separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter , Value )


PAM.piazza.400.2.1 <- bind_rows (Inv.param.400.2, LB.param.400.2 )

PAM.piazza.400.2.1  <- PAM.piazza.400.2.1 %>%
                       dplyr::mutate (ppfd.act = 200)


PAM.piazza.400.2.1 <- PAM.piazza.400.2.1 %>%
  dplyr::rename (sampling = ambiente) %>%
  dplyr::mutate (protocol = "piazza.400.2" )


write_delim (PAM.piazza.400.2.1, path ="./Data/procdata/PAM.piazza.400.2.1.txt" , 
             delim = ",", na = "NA")

PAM.param.pza800 <- read_delim (file ="./Data/procdata/pre.PAM.param.pza800.1.txt", 
                                delim =",", quote = "\"",
                                escape_backslash = FALSE,
                                escape_double = TRUE, 
                                col_names = TRUE, 
                                col_types = NULL, 
                                locale = default_locale(), 
                                na = "NA")


PAM.param.pza800.1 <- PAM.param.pza800 %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") 

head (PAM.param.pza800.1)

Inv.param.800 <- PAM.param.pza800.1 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "800" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter ,  Value )


LB.param.800  <-  PAM.param.pza800.1  %>%
  dplyr::filter (ambiente != "Inv")

LB.param.800  <- LB.param.800  %>% 
  tidyr::separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter , Value )


PAM.piazza.800.1 <- bind_rows (Inv.param.800, LB.param.800 )

PAM.piazza.800.1  <- PAM.piazza.800.1 %>%
                     dplyr::mutate (ppfd.act = 1.7e3)


PAM.piazza.800.1 <- PAM.piazza.800.1 %>%
  dplyr::rename (sampling = ambiente) %>%
  dplyr::mutate (protocol = "piazza.800" )


write_delim (PAM.piazza.800.1, path ="./Data/procdata/PAM.piazza.800.1.txt" , 
             delim = ",", na = "NA")

#### Piazza 800.2

PAM.param.pza800.2 <- read_delim (file ="./Data/procdata/pre.PAM.param.pza800.2.1.txt", 
                                  delim =",", quote = "\"",
                                  escape_backslash = FALSE,
                                  escape_double = TRUE, 
                                  col_names = TRUE, 
                                  col_types = NULL, 
                                  locale = default_locale(), 
                                  na = "NA")


PAM.param.pza800.2.1 <- PAM.param.pza800.2 %>%
  dplyr::mutate (x = pot)%>%
  tidyr::separate (x , c("ambiente","ppfd","genotype", "rep", NA), sep ="_") 

head (PAM.param.pza800.2.1)

Inv.param.800.2 <- PAM.param.pza800.2.1 %>%
  dplyr::filter ( ambiente == "Inv") %>%
  dplyr::mutate ( ppfd = "800.2" ) %>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter ,  Value )


LB.param.800.2  <-  PAM.param.pza800.2.1  %>%
  dplyr::filter (ambiente != "Inv")

LB.param.800.2  <- LB.param.800.2  %>% 
  tidyr::separate(ppfd, c("ppfd", NA))%>%
  dplyr::select (ambiente,  ppfd, genotype, rep, pot, Parameter , Value )


PAM.piazza.800.2.1 <- bind_rows (Inv.param.800.2, LB.param.800.2 )

PAM.piazza.800.2.1  <- PAM.piazza.800.2.1 %>%
                       dplyr::mutate (ppfd.act = 850)

PAM.piazza.800.2.1 <- PAM.piazza.800.2.1 %>%
  dplyr::rename (sampling = ambiente) %>%
  dplyr::mutate (protocol = "piazza.800.2" )


write_delim (PAM.piazza.800.2.1, path ="./Data/procdata/PAM.piazza.800.2.1.txt" , 
             delim = ",", na = "NA")


################################ Aca reingreso ### 
PAM.param.pza400 <- read_delim (file ="./Data/procdata/PAM.piazza.400.1.txt", 
                                  delim =",", quote = "\"",
                                  escape_backslash = FALSE,
                                  escape_double = TRUE, 
                                  col_names = TRUE, 
                                  col_types = NULL, 
                                  locale = default_locale(), 
                                  na = "NA")


PAM.param.pza400.2 <- read_delim (file ="./Data/procdata/PAM.piazza.400.2.1.txt", 
                                  delim =",", quote = "\"",
                                  escape_backslash = FALSE,
                                  escape_double = TRUE, 
                                  col_names = TRUE, 
                                  col_types = NULL, 
                                  locale = default_locale(), 
                                  na = "NA")


PAM.param.pza800 <- read_delim (file ="./Data/procdata/PAM.piazza.800.1.txt", 
                                delim =",", quote = "\"",
                                escape_backslash = FALSE,
                                escape_double = TRUE, 
                                col_names = TRUE, 
                                col_types = NULL, 
                                locale = default_locale(), 
                                na = "NA")


PAM.param.pza800.2 <- read_delim (file ="./Data/procdata/PAM.piazza.800.2.1.txt", 
                                  delim =",", quote = "\"",
                                  escape_backslash = FALSE,
                                  escape_double = TRUE, 
                                  col_names = TRUE, 
                                  col_types = NULL, 
                                  locale = default_locale(), 
                                  na = "NA")

 # 


#dt.curve = PAM.curve.GT.1.s3
#dt.param = PAM.param.GT.1.s3.a
#time.quenching = 	435.05
#sampling ="GT.1"
#protocol = "s3"
#rep=TRUE


#dt = PAM.param.pza400
#sampling = NULL
#protocol = "piazza.400"
#ppfd.act = 403 



#simondi_analisis_Fo


############# Esta es la funcion
run_simondi_analisis_Fo_pzza <-  function ( dt = NULL, sampling = NULL, protocol = NULL, ppfd.act = NULL) {
  
  dir.create (file.path ("Figures", "Plots.phi.simondi.Fo"), showWarnings = FALSE)
  
  prt <- protocol
  param.protocol <- dt %>%
                    dplyr::filter (protocol == prt)
  
  list.pot <- unique (param.protocol$pot)
  
  dt.pot <- lapply(list.pot, function (filt.pot){
    
    #filt.pot = "LB_800.8_HHK2C5_1_Param.asc"
    
    pot.par <- dt %>%
               dplyr::filter (pot == filt.pot )
    
    pot.par <- pot.par %>%
               dplyr::mutate (ord = as.numeric(row.names (pot.par))) %>%
               dplyr::select (ord, everything())
    
    id <- pot.par %>%
          dplyr::select (c(sampling, pot, protocol, genotype, rep , ppfd.act )) %>%
          dplyr::select (sampling, protocol,ppfd.act, genotype,rep,  pot  )
    
    idx <- str_c (unique (pot.par$pot),
                  unique (pot.par$ppfd.act), 
                  sep="_")
    
    ## parametros del PAM
    fluor.min <- pot.par %>%
                 dplyr::filter (ord < 3)  %>%
                 dplyr::summarise (Fo = min (Value))
    
    fluor.max <- pot.par %>%
                dplyr::filter (ord == 3) %>%
                dplyr::mutate (Fm = Value)%>%
                dplyr::select (Fm)
    
    fluor.time <- pot.par %>%
                  dplyr::filter (Parameter == "Fs") 
    

    fluor.max.prim <- pot.par %>%
                      dplyr::filter (Parameter == "Fm'")

    
    fluor.min.prim <- pot.par %>%
                      dplyr::filter (Parameter == "Fo'") %>%
                      dplyr::filter (ord <= 39)
    
    
    ##### calculamos parametros Fo variables ####
    
    phi.ps2 <- tibble (phiPS2 = (fluor.max.prim$Value - fluor.time$Value)/fluor.max.prim$Value)
    
    Q.p.fo <- tibble ( qp.fo = (fluor.max.prim$Value - fluor.time$Value)/(fluor.max.prim$Value - fluor.min.prim$Value ))
    S.p.fo <- tibble ( sp.fo = (fluor.max.prim$Value - fluor.min.prim$Value)/fluor.max.prim$Value)
    
    
    
    phi.NPQ <- tibble (phiNPQ = (fluor.time$Value/fluor.max.prim$Value) * ((fluor.max$Fm - fluor.max.prim$Value)/ fluor.max$Fm ))
    
    NPQ.1.fo   <- tibble (npq1.fo = ((fluor.time$Value - fluor.min.prim$Value)/fluor.max.prim$Value) * ((fluor.max$Fm - fluor.max.prim$Value)/fluor.max$Fm ))
    
    NPQ.2.fo   <- tibble (npq2.fo = (fluor.min.prim$Value/fluor.max.prim$Value) *((fluor.max$Fm - fluor.max.prim$Value)/fluor.max$Fm))
    
    
    phi.NO <- tibble (phiNO =fluor.time$Value/fluor.max$Fm)
    
    NO.1.fo <- tibble (no1.fo = (fluor.time$Value - fluor.min.prim$Value) /fluor.max$Fm )
    ### esto es centro de reaccion
    
    
    NO.2.fo <- tibble (no2.fo = fluor.min.prim$Value/fluor.max$Fm)
    ### esto es complejo antena
    
    sum.phi.fo <-  tibble (sphi = phi.ps2$phiPS2 + phi.NPQ$phiNPQ + phi.NO$phiNO)
    
    quenching.analisis.fo <- bind_cols(sum.phi.fo, phi.ps2, Q.p.fo, S.p.fo, phi.NPQ, NPQ.1.fo, NPQ.2.fo, phi.NO, NO.1.fo, NO.2.fo)
    
    quenching.analisis.fo <- quenching.analisis.fo %>%
      dplyr::mutate  (ciclo = row.names (quenching.analisis.fo)) 
    
    x <- length(unique (quenching.analisis.fo$ciclo))-1
    x1 <- c(1:x * 30) + 60
    time <- c(60,x1 )
    
    quenching.analisis.fo <- quenching.analisis.fo %>%
      dplyr::mutate (time.min = time/60)
    
    
    quenching.analisis.fo.1 <- quenching.analisis.fo %>%
      dplyr::select (-c (sphi, qp.fo, sp.fo ))
    
    
    quenching.analisis.fo.1t <- quenching.analisis.fo.1 %>%
      pivot_longer(-c(ciclo, time.min), names_to = "quantum.yield" ,values_to = "phi")
    
    part.simondi.fo <- c("phiPS2","npq1.fo","npq2.fo","no1.fo","no2.fo")
    
    quenching.analisis.fo.1ts <- quenching.analisis.fo.1t  %>%
      dplyr::filter (quantum.yield %in%  part.simondi.fo)%>%
      dplyr::mutate (quantum.yield = factor(quantum.yield, levels = c("phiPS2","npq1.fo", "npq2.fo","no1.fo", "no2.fo")))%>%
      dplyr::mutate (phi = round (phi,2))
    
    phi.plot.simondi.fo <- ggbarplot(quenching.analisis.fo.1ts , "time.min", "phi",title =idx,
                                     ylab="time (min)",
                                     fill = "quantum.yield", color = "quantum.yield", 
                                     palette = c("palevioletred3","steelblue3","steelblue4" ,"wheat3", "wheat4" ),
                                     label = TRUE, lab.col = "white", lab.pos = "in")
    
    
    print (phi.plot.simondi.fo)
    
    ggexport (phi.plot.simondi.fo, filename = str_c("./Figures/Plots.phi.simondi.Fo/",prt,"_",sampling,"_" , "phi.plot.Fo", "_",idx,".tiff"),
              width = 700, height = 500)
    
  
    xx <- id [1,]
    quenching.analisis.2.fo <- cbind (xx, quenching.analisis.fo)
    
    
    return (quenching.analisis.2.fo)
    
  })# aca termina dt.pot
  
  df.ss <- bind_rows (dt.pot )
  return (df.ss)
} ### run_simondi_analisis_GT 
    
    
    
  

#######################
unique (PAM.param.pza400$protocol )
unique (PAM.param.pza400.2$protocol )

phi.PAM.param.pza400.Fo <- run_simondi_analisis_Fo_pzza (dt = PAM.param.pza400, sampling = NULL, protocol = "piazza.400", ppfd.act = 403 )


phi.PAM.param.pza400.2.Fo <- run_simondi_analisis_Fo_pzza (dt = PAM.param.pza400.2, sampling = NULL, protocol = "piazza.400.2", ppfd.act = 200 )


phi.PAM.param.pza800.Fo <- run_simondi_analisis_Fo_pzza (dt = PAM.param.pza800, sampling = NULL, protocol = "piazza.800", ppfd.act = 1700 )


phi.PAM.param.pza800.2.Fo <- run_simondi_analisis_Fo_pzza (dt = PAM.param.pza800.2, sampling = NULL, protocol = "piazza.800.2", ppfd.act = 850 )


phi.PAM.piazzas.Fo <- bind_rows (phi.PAM.param.pza400.Fo, phi.PAM.param.pza400.2.Fo , phi.PAM.param.pza800.Fo , phi.PAM.param.pza800.2.Fo )



write_delim (phi.PAM.piazzas.Fo, file ="./Data/procdata/phi.PAM.piazzas.Fo.txt" , 
             delim = ",", na = "NA")

############## hasta el 20/01/2020 ########################################


Inv.phi.400 <- run_quenching_analisis_pzza (dt = PAM.piazza.400.1, 
                                             sampling = "Inv",
                                             protocol = "piazza.400", 
                                             ppfd.act = 403)



LB.PAM.piazza.400.1 <- PAM.piazza.400.1 %>%
                       dplyr::filter ( sampling != "Inv" )


LB.phi.400 <- run_quenching_analisis_pzza (dt = PAM.piazza.400.1, 
                                            sampling = "Inv",
                                            protocol = "piazza.400", 
                                            ppfd.act = 403)





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
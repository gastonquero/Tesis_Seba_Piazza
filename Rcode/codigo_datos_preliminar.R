##############################################################################
# Codigo para la funcion de PAM
# 
# Gaston Quero - Sebastian Simondi                                
# 14/09/2020
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

files.curva <- dir(str_c("./Data/rawdata/preliminares/Resultados/curvas"), pattern = "*.asc")

PAM.curves <- bind_rows (lapply(files.curva, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/preliminares/Resultados/curvas/", filt.raw) , 
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
  

X <- lapply(files.curva, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves %>%
         dplyr::filter (pot == filt.pot)
         
  run.plot.curve (dt=pot, sampling = "preliminar" )    
  
  
})


PAM.curves.preliminar  <- PAM.curves %>%
                          dplyr::mutate (x = pot)%>%
                          tidyr::separate(x , c("protocol", "genotype", "rep", NA), sep ="_") %>%
                          dplyr::select (-X3)

write_delim (PAM.curves.preliminar, path ="./Data/procdata/curvas_preliminares.txt" , 
             delim = ",", na = "NA")


## Aca cargo los datos de los parametros
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
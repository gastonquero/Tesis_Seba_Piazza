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
#C:\Users\Usuario\OneDrive\Documentos\Tesis_Seba_Piazza\Data\rawdata\Piazza_400\parametros_pza400
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

### aca reingreso los datos 
PAM.param.pza400 <- read_delim (file ="./Data/procdata/pre.PAM.param.pza400.1.txt", 
                         delim =",", quote = "\"",
                         escape_backslash = FALSE,
                         escape_double = TRUE, 
                         col_names = TRUE, 
                         col_types = NULL, 
                         locale = default_locale(), 
                         na = "NA")



unique (PAM.param.pza400$pot)


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




######
### hay que ir por protocolo
sampling = "preliminar_"

#param.s3 <- PAM.param.1 %>%
            #dplyr::filter ( protocol == "Simondi3")

list.pot.pza400 <- unique (PAM.piazza.400.1$pot)
unique (PAM.param.1$protocol)

dt = PAM.piazza.400.1
#head ( PAM.piazza.400.1)
protocol = "piazza.400"
sampling = "Inv"
ppfd.act = 403
  
run_quenching_analisis_pzza <-  function ( dt = NULL, sampling = NULL, protocol = NULL, ppfd.act = NULL) {
  
  prt <- protocol
  param.protocol <- dt %>%
                    dplyr::filter (protocol == prt)
  
  list.pot <- unique (param.protocol$pot)
  
  df.pot <- bind_rows (lapply(list.pot, function (filt.pot){
    
    #filt.pot = "Inv_1_HHK2C5_4_Param.asc"
    
    pot.par <- dt %>%
               dplyr::filter (pot == filt.pot )
    
    pot.par <- pot.par %>%
               dplyr::mutate (ord = as.numeric(row.names (pot.par))) %>%
               dplyr::select (ord, everything())
    
    id <- pot.par %>%
          dplyr::select (c(pot, protocol, genotype, rep , ppfd.act )) %>%
          dplyr::select (protocol,ppfd.act, genotype,rep,  pot  )
    
    ## parametros inciales
    fluor.min <- pot.par %>%
                 dplyr::filter (ord < 3) %>%
                 dplyr::summarise (Fo = mean (Value))
    
    fluor.max <- pot.par %>%
                dplyr::filter (ord == 3) %>%
                dplyr::mutate (Fm = Value)%>%
                dplyr::select (Fm)
    
    fluor.max.prim <- pot.par %>%
                      dplyr::filter (Parameter == "Fm'" )

    
    fluor <- pot.par  %>%
             dplyr::filter (Parameter == "Fs") 
    
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
    
    phi.plot <- ggbarplot(quenching.analisis.1 , "time.min", "phi",title =str_c("ppfd.act.",ppfd.act,unique (id$pot), sep="_"),
                          ylab="time (min)",
                          fill = "quantum.yield", color = "quantum.yield", 
                          palette = c("gray18","gray38" ,"gray58"),
                          label = TRUE, lab.col = "white", lab.pos = "in")
    
    ggexport (phi.plot, filename = str_c("./Figures/Plots.PAM/",prt,"phi.plot", unique(id$pot), ".tiff"),
              width = 700, height = 500)
    
    print (phi.plot)
    xid <-id [1,] 
    
    quenching.analisis.2 <- cbind (quenching.analisis, xid)
    
    return (quenching.analisis.2 )
  }))
  
  return (df.pot)
}

Inv.PAM.piazza.400.1 <- PAM.piazza.400.1 %>%
                        dplyr::filter ( sampling == "Inv" )

  
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
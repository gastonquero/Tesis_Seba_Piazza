######################################################################################
# Codigo para el analisis de los Datos de la Tesis de Sebastian Piazza               #
#                                                                                    #
# Gaston Quero - Sebastian Simondi                                                   #          
# 8/7/2021                                                                          #
######################################################################################

#####
# los datos correspoden a .......

getwd ()
setwd ( "R:/Tesis_Seba_Piazza" )

#setwd("C:/Users/Usuario/OneDrive/Documentos/Paper_fotosintesis_dinámica")

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
library("Hmisc")
library("PerformanceAnalytics")
library (GGally)
library(ggcorrplot)


### seteo del directorio donde estan los datos

dir.data.base <- ("C:/Users/Usuario/OneDrive/Documentos/Paper_fotosintesis_dinámica")

## SP.1
Phi.PSII.SP.1 <- read_delim (file =str_c(dir.data.base,"/Data/procdata/SP_1/Phi.PSII.SP.1.txt"), 
                             delim =";", quote = "\"", escape_backslash = FALSE,
                             escape_double = TRUE, col_names = TRUE, col_types = NULL,
                             na = "NA")


Phi.PSII.SP.1a <- Phi.PSII.SP.1 %>%
                  tidyr::separate (id , c("ensayo","ambiente","genotipo","cond.hidr", "pos", "hoja","protocolo", "date"), sep="_") 




### unifico las matrices 

unique (Phi.PSII.SP.1a$ensayo)
unique (Phi.PSII.SP.1a$ambiente)
unique (Phi.PSII.SP.1a$genotipo)
unique (Phi.PSII.SP.1a$pos)
unique (Phi.PSII.SP.1a$protocolo)


# protocolos cambio a actinicas
Phi.PSII.SP.1a <- Phi.PSII.SP.1a  %>%
                  dplyr::mutate (ppfd.act = protocolo)
                       

Phi.PSII.SP.1a$ppfd.act [Phi.PSII.SP.1a$ppfd.act == "pz400.2"]      <- 200
Phi.PSII.SP.1a$ppfd.act [Phi.PSII.SP.1a$ppfd.act == "pz400"]        <- 400
Phi.PSII.SP.1a$ppfd.act [Phi.PSII.SP.1a$ppfd.act == "pz800.2"]      <- 850
Phi.PSII.SP.1a$ppfd.act [Phi.PSII.SP.1a$ppfd.act == "pz800"]        <- 1700



SP1.LB400  <- Phi.PSII.SP.1a %>%
              dplyr::filter (str_detect(filt.pot, "LB_400")) %>%
              dplyr::mutate (ppfd.amb = 400)

SP1.LB800  <- Phi.PSII.SP.1a %>%
              dplyr::filter (str_detect(filt.pot, "LB_800")) %>%
              dplyr::mutate (ppfd.amb = 800)

SP1.Inv  <- Phi.PSII.SP.1a %>%
            dplyr::filter (str_detect(filt.pot, "Inv")) %>%
            dplyr::mutate (ppfd.amb = 450)



Phi.PSII.SP.1b <- bind_rows (SP1.LB400 , SP1.LB800 ,SP1.Inv) 


Phi.PSII.SP.1b  <- Phi.PSII.SP.1b  %>%
                   dplyr::select (ensayo, ambiente,ppfd.amb,protocolo, ppfd.act, genotipo ,cond.hidr, pos,  hoja, everything() )

write_excel_csv (Phi.PSII.SP.1b , file ="./Data/procdata/Phi.PSII.SP.1b.csv",
                  na = "NA",
                  append = FALSE,
                  delim = ";",
                  quote_escape = "double",
                  eol = "\n" )

write_delim (Phi.PSII.SP.1b , file ="./Data/procdata/Phi.PSII.SP.1b.txt", 
             delim = ";", na = "NA")


# saco el ultimo ciclo #### 

pz400.quenching <- Phi.PSII.SP.1b %>%
                   dplyr::filter (protocolo == "pz400" ) %>%
                   dplyr::filter (ciclo == max (ciclo) - 1)


pz400.2.quenching <- Phi.PSII.SP.1b %>%
                     dplyr::filter (protocolo == "pz400.2" ) %>%
                     dplyr::filter (ciclo == max (ciclo) - 1)

pz800.quenching <- Phi.PSII.SP.1b %>%
                   dplyr::filter (protocolo == "pz800" ) %>%
                   dplyr::filter (ciclo == max (ciclo) - 1)


pz800.2.quenching <- Phi.PSII.SP.1b %>%
                     dplyr::filter (protocolo ==  "pz800.2" ) %>%
                     dplyr::filter (ciclo == max (ciclo) - 1)


Phi.PSII.SP.1.ciclo.q <- bind_rows ( pz400.quenching,
                                          pz400.2.quenching,
                                          pz800.quenching,
                                          pz800.2.quenching)


write_excel_csv (Phi.PSII.SP.1.ciclo.q, file ="./Data/procdata/Phi.PSII.SP.1.ciclo.q.csv",
                 na = "NA",
                 append = FALSE,
                 delim = ";",
                 quote_escape = "double",
                 eol = "\n" )

write_delim (Phi.PSII.SP.1.ciclo.q, file ="./Data/procdata/Phi.PSII.SP.1.ciclo.q.txt", 
             delim = ";", na = "NA")

##################################### Plot ###
# filtro los datos fuera de rango

phiNPQslow.neg <- Phi.PSII.SP.1.ciclo.q %>%
                  dplyr::filter (phiNPQslow < 0 ) %>%
                  dplyr::arrange (phiNPQslow )

Phi.PSII.SP.1.ciclo.q.1  <- Phi.PSII.SP.1.ciclo.q %>%
                            dplyr::filter (phiNPQslow > -0.009 ) %>%
                            dplyr::filter ( qs.fo < 1 )


Phi.PSII.SP.1.ciclo.q.1$phiNPQslow[Phi.PSII.SP.1.ciclo.q.1$phiNPQslow < 0] <- 0


Phi.PSII.SP.1.ciclo.q.1 <- Phi.PSII.SP.1.ciclo.q.1 %>%
                           dplyr::mutate (trat = str_c (ambiente, protocolo, genotipo)) %>%
                           dplyr::mutate (amb.ppfd = str_c (ambiente, ppfd.amb)) 
  
##### 
unique ( Phi.PSII.SP.1.ciclo.q.1$amb.ppfd  )

LB800.1700 <- Phi.PSII.SP.1.ciclo.q.1 %>%
              dplyr::filter (amb.ppfd ==  "LB800") %>%
              dplyr::filter ( ppfd.act  == 1700)

potX <- unique (LB800.1700$filt.pot)

Phi.PSII.SP.1.ciclo.q.1 <-  Phi.PSII.SP.1.ciclo.q.1 %>%
                             dplyr::filter (!filt.pot %in% potX)




# garfico de correlaciones 

phis <- c ("phiPS2", "qp.fo" , "qs.fo" , "phiNPQ", "phiNO", "phiNO.psII" ,
            "phiNO.basal" , "phiNPQfast", "phiNPQslow" ) 


phis.cor <- Phi.PSII.SP.1.ciclo.q.1 %>%
            dplyr::select (all_of(phis))

ggpairs (phis.cor)

corr <- round (cor (phis.cor), 2)

p.mat <- cor_pmat (phis.cor )

ggcorrplot(corr, hc.order = TRUE, type = "lower",
           colors = c("darkred", "white", "blue"),insig ="blank",
           lab = TRUE,p.mat = p.mat, sig.level =0.05)

#### corremos el modelo 
## model 
lm.phiPS2 <- lm (phiPS2 ~  ambiente + protocolo + genotipo + 
                       ambiente * protocolo * genotipo 
                     , data= Phi.PSII.SP.1.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiPS2$fitted.values, lm.phiPS2$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiPS2")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiPS2$residuals, main="phiPS2",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiPS2$residuals),
              sd (lm.phiPS2$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiPS2$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiPS2")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiPS2$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiPS2 ~ trat, data = Phi.PSII.SP.1.ciclo.q.1)

# ANAVA
(anava.lm.phiPS2  <- anova (lm.phiPS2))

############## qp ######
lm.qp.fo <- lm (qp.fo ~  ambiente + protocolo + genotipo + 
                   ambiente * protocolo * genotipo 
                 , data= Phi.PSII.SP.1.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.qp.fo$fitted.values, lm.qp.fo$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="qp.fo")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.qp.fo$residuals, main="qp.fo",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.qp.fo$residuals),
              sd (lm.qp.fo$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.qp.fo$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="qp.fo")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.qp.fo$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( qp.fo ~ trat, data = Phi.PSII.SP.1.ciclo.q.1)

# ANAVA
(anava.lm.qp.fo  <- anova (lm.qp.fo))


####
Phi.PSII.SP.1.ciclo.q.1$protocolo  <- factor(Phi.PSII.SP.1.ciclo.q.1$protocolo, 
                                      levels = c("pz400.2", 
                                               "pz400",
                                              "pz800.2", "pz800"))

Phi.PSII.SP.1.ciclo.q.1$ppfd.act  <- factor(Phi.PSII.SP.1.ciclo.q.1$ppfd.act, 
                                             levels = c("200", 
                                                        "400",
                                                        "850", "1700"))
### hacer los sactter

ggscatter (  Phi.PSII.SP.1.ciclo.q.1, x = "qp.fo", y = "phiPS2",
             ylim=c(0,1),
             #xlim=c(0,1),
             facet.by =  "ppfd.act",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "amb.ppfd",
             size =5,
             palette = c("darkorange","darkred", "navyblue")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)




############## NPQ ######
pot_X <-  Phi.PSII.SP.1.ciclo.q.1 %>%
          dplyr::filter (filt.pot == "pz400_LB_400.8_HHK1C3_4_Param.asc")





#Phi.PSII.SP.1.ciclo.q.1 <- Phi.PSII.SP.1.ciclo.q.1  %>%
 #                          dplyr::filter (filt.pot != "pz400_LB_400.8_HHK1C3_4_Param.asc" )



lm.phiNPQ <- lm (1/log(phiNPQ) ~  ambiente + protocolo + genotipo + 
                  ambiente * protocolo * genotipo 
                , data= Phi.PSII.SP.1.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiNPQ$fitted.values, lm.phiNPQ$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiNPQ")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiNPQ$residuals, main="phiNPQ",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiNPQ$residuals),
              sd (lm.phiNPQ$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiNPQ$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiNPQ")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiNPQ$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( 1/log(phiNPQ) ~ trat, data = Phi.PSII.SP.1.ciclo.q.1)

# ANAVA
(anava.lm.phiNPQ  <- anova (lm.phiNPQ))


############## phiNPQ fast ##################33
head ( Phi.PSII.SP.1.ciclo.q.1 )

############## NPQ ######
pot_X <-  Phi.PSII.SP.1.ciclo.q.1 %>%
  dplyr::filter (filt.pot == "pz400_LB_400.8_HHK1C3_4_Param.asc")





#Phi.PSII.SP.1.ciclo.q.1 <- Phi.PSII.SP.1.ciclo.q.1  %>%
#                          dplyr::filter (filt.pot != "pz400_LB_400.8_HHK1C3_4_Param.asc" )



lm.phiNPQfast <- lm ( 1/log(phiNPQfast) ~  ambiente + protocolo + genotipo + 
                   ambiente * protocolo * genotipo 
                 , data= Phi.PSII.SP.1.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiNPQfast$fitted.values, lm.phiNPQfast$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiNPQ")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiNPQfast$residuals, main="phiNPQ",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiNPQfast$residuals),
              sd (lm.phiNPQfast$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiNPQfast$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiNPQ")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiNPQfast$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest (1/log(phiNPQfast) ~ trat, data = Phi.PSII.SP.1.ciclo.q.1)

# ANAVA
(anava.lm.phiNPQfast  <- anova (lm.phiNPQfast))

  
  
ggscatter ( Phi.PSII.SP.1.ciclo.q.1, x = "phiNPQfast", y = "phiNPQ",
             ylim=c(0,1),
             #xlim=c(0,1),
             facet.by = "ppfd.act",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "amb.ppfd",
             size =5,
             palette = c("darkorange","darkred", "navyblue")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)

head (Phi.PSII.SP.1.ciclo.q.1)

ggscatter ( Phi.PSII.SP.1.ciclo.q.1, x = "qp.fo", y = "phiNPQfast",
            ylim=c(0,1),
            #xlim=c(0,1),
            facet.by = "ppfd.act",
            #add="reg.line",
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            color= "amb.ppfd",
            size =5,
            palette = c("darkorange","darkred", "navyblue")
            
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)



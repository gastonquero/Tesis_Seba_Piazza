########################################################
#  Analisis de Datos Luz actinica PAM                  #
#                                                      #
# Gaston Quero - Sebastian Fernandez                   #
#           6-03-2017                                  #
########################################################

getwd()
setwd("C:/Users/Usuario/OneDrive/Documentos/Paper_fotosintesis_dinámica")
# seteo el directorio  #### Aca se cambia segun donde tengas los datos 
                          # Esta linea se corre una sola vez en cada caso, y se 
                          # desabilita la linea del directorio de otro usuario
                          # esto lo usamos asi hasta que trabajemos en Git



# Paquetes 
library(lme4)
library(lmerTest)
library(nlme)
library(car)
library("ggplot2")       
library("lattice")
library("latticeExtra")
library(multcompView)
library(dplyr)
library(plyr)
library(xtable)
library(tidyverse)
library (emmeans)
library("qtl")
library(stringr)
library(data.table)
library(svMisc)
library(ggpubr)
library("ggsci")
library("FactoMineR")
library("factoextra")
library("corrplot")
library(Matrix)

######### distribucion espectral  ####################################
#########  Ambiente Azul y Rojo  #####################################
#### Actinica 25 lampara nueva  ########################################################

# directorio sebfer
#[dist_espectral]= uW/(cm2.nm)
#dist_espectral_uw_raw <- read.table("./lampara_hpit_AbsoluteIrradiance_0002.txt", header = FALSE, sep = "\t",dec = ",", skip=14)

# Diretorio GQ


## GT.1 
dist_espectral_uw_GT.1a_raw <- read.table("./Data/rawdata/ambientes_lumínicos_cebada/cebada.invernaculo_AbsoluteIrradiance_12-57-51-130.txt",
                                         header = FALSE, sep = "\t",dec = ",", skip=14)

colnames(dist_espectral_uw_GT.1a_raw) <- c("lambda", "uW_nm.cm2")

dist_espectral_uw_GT.1a_raw  <- dist_espectral_uw_GT.1a_raw%>%
                                dplyr::mutate (punto= "X")



write_delim (dist_espectral_uw_GT.1a_raw, path= "./Data/procdata/dist_espectral_uw_GT.1_raw.txt", delim = ",", na = "NA")


## GT.2 
dist_espectral_uw_GT.2_raw <- read.table("./Data/rawdata/ambientes_lumínicos_cebada/GT.2_AbsoluteIrradiance_00001.txt",
                                          header = FALSE, sep = "\t",dec = ",", skip=14)

#View (dist_espectral_uw_raw)
colnames(dist_espectral_uw_GT.2_raw) <- c("lambda", "uW_nm.cm2")

dist_espectral_uw_GT.2_raw  <- dist_espectral_uw_GT.2_raw  %>%
                               dplyr::mutate (punto = "X")


write_delim (dist_espectral_uw_GT.2_raw, path= "./Data/procdata/dist_espectral_uw_GT.2_raw.txt", delim = ",", na = "NA")






## GT.3 
dist_espectral_uw_GT.3_raw <- read.table("./Data/rawdata/ambientes_lumínicos_cebada/GT.3_AbsoluteIrradiance_00001.txt",
                                         header = FALSE, sep = "\t",dec = ",", skip=14)

#View (dist_espectral_uw_raw)
colnames(dist_espectral_uw_GT.3_raw) <- c("lambda", "uW_nm.cm2")

dist_espectral_uw_GT.3_raw  <- dist_espectral_uw_GT.3_raw  %>%
                               dplyr::mutate (punto = "X")


write_delim (dist_espectral_uw_GT.3_raw, path= "./Data/procdata/dist_espectral_uw_GT.3_raw.txt", delim = ",", na = "NA")



## GT.4 
dist_espectral_uw_GT.4_raw <- read.table("./Data/rawdata/ambientes_lumínicos_cebada/GT.4_AbsoluteIrradiance_01.txt",
                                         header = FALSE, sep = "\t",dec = ",", skip=14)

#View (dist_espectral_uw_raw)
colnames(dist_espectral_uw_GT.4_raw) <- c("lambda", "uW_nm.cm2")

dist_espectral_uw_GT.4_raw  <- dist_espectral_uw_GT.4_raw  %>%
   dplyr::mutate (punto = "X")

write_delim (dist_espectral_uw_GT.4_raw, path= "./Data/procdata/dist_espectral_uw_GT.4_raw.txt", delim = ",", na = "NA")




# Construccion de la matriz de lambdas
d1 <- tibble ( bandwidth = "D1", l1 = 400, l2= 425)
d2 <- tibble ( bandwidth = "D2", l1 = 425, l2= 490)
d3 <- tibble ( bandwidth = "D3", l1 = 490, l2= 560)
d4 <- tibble ( bandwidth = "D4", l1 = 560, l2= 585)
d5 <- tibble ( bandwidth = "D5", l1 = 585, l2= 640)
d6 <- tibble ( bandwidth = "D6", l1 = 640, l2= 700)
dGT2 <- tibble ( bandwidth = "GT2", l1 = 420, l2= 649)
dGT3 <- tibble ( bandwidth = "GT3", l1 = 415, l2= 584)
dGT4 <- tibble ( bandwidth = "GT4", l1 = 550, l2= 700)
dPAR <- tibble ( bandwidth = "PAR", l1 = 400, l2= 700)
dintegr <- tibble ( bandwidth = "total", l1 = 380, l2= 780)

lambdas <- bind_rows (d1, d2, d3, d4, d5, d6, dGT2, dGT3, dGT4,
                      dPAR,dintegr )


####
amb.luminico.GT.1 <- run_spectrum (dt = dist_espectral_uw_GT.1a_raw, lambdas = lambdas, id.survey="GT.1", px="X" )

amb.luminico.GT.1  <- amb.luminico.GT.1 %>%
                     dplyr::mutate (relevamiento = "amb.luminico.GT.1")


amb.luminico.GT.2 <- run_spectrum (dt = dist_espectral_uw_GT.2_raw, lambdas = lambdas, id.survey="GT.2", px="X" )

amb.luminico.GT.2  <- amb.luminico.GT.2 %>%
                      dplyr::mutate (relevamiento = "amb.luminico.GT.2")


amb.luminico.GT.3 <- run_spectrum (dt = dist_espectral_uw_GT.3_raw, lambdas = lambdas, id.survey="GT.3", px="X" )

amb.luminico.GT.3  <- amb.luminico.GT.3 %>%
                      dplyr::mutate (relevamiento = "amb.luminico.GT.3")

amb.luminico.GT.4 <- run_spectrum (dt = dist_espectral_uw_GT.4_raw, lambdas = lambdas, id.survey="GT.4", px="X" )

amb.luminico.GT.4  <- amb.luminico.GT.4 %>%
                      dplyr::mutate (relevamiento = "amb.luminico.GT.4")




calidades.luminicas.fot.cebada <- bind_rows(amb.luminico.GT.1,
                                            amb.luminico.GT.2,
                                            amb.luminico.GT.3,
                                            amb.luminico.GT.4)


write_delim (calidades.luminicas.fot.cebada, path= "./Data/procdata/calidades.luminicas.fot.cebada.txt", delim = ",", na = "NA")

#### codigo para hacer la figura de los ambientes

amb_GT.1  <- dist_espectral_uw_GT.1a_raw %>%
             dplyr::mutate (amb = "GT.1")%>%
             dplyr::select (-punto)

amb_GT.2  <- dist_espectral_uw_GT.2_raw %>%
             dplyr::mutate (amb = "GT.2")%>%
             dplyr::select (-punto)
             
amb_GT.3  <- dist_espectral_uw_GT.3_raw %>%
             dplyr::mutate (amb = "GT.3")%>%
             dplyr::select (-punto)        
             
amb_GT.4  <- dist_espectral_uw_GT.4_raw %>%
             dplyr::mutate (amb = "GT.4")%>%
             dplyr::select (-punto)             
             

amb_GT <- bind_rows (amb_GT.1,amb_GT.2,amb_GT.3,amb_GT.4)             
             

plot (uW_nm.cm2 ~ lambda, 
      ylim = c(0,150),
      xlim=  c(400, 750),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c ( "power",unique(amb_GT.1$amb) ,sep="_",collapse = TRUE),
      data= amb_GT.1) 

abline ( v=420, lty=2, col="gray48" )
abline ( v=649 , lty=2, col="gray48" )


x1.pwr <- 619
x2.pwr <- 1279


amb_GT.1$lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]

with (amb_GT.1, polygon(x=c(lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]),
                       y= c(0, uW_nm.cm2[x1.pwr:x2.pwr], 0),border ='gray28',
                       col=scales::alpha( 'gray28',.3)))

text(x = 500 , 
     y= 100 ,
     label = "intervalo = 420-649 nm;
              PPFD = 349 um/m2/s;
              power= 78 W/m2",
     pos=3,
     cex = 0.9,  col="black")



###

plot (uW_nm.cm2 ~ lambda, 
      ylim = c(0,150),
      xlim=  c(400, 750),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c ( "power",unique(amb_GT.2$amb) ,sep="_",collapse = TRUE),
      data= amb_GT.2) 

abline ( v=420, lty=2, col="gray48" )
abline ( v=649 , lty=2, col="gray48" )


x1.pwr <- 619
x2.pwr <- 1279


#amb_GT.1$lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]

with (amb_GT.2, polygon(x=c(lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]),
                        y= c(0, uW_nm.cm2[x1.pwr:x2.pwr], 0),border ='gray28',
                        col=scales::alpha( 'gray28',.3)))

text(x = 500 , 
     y= 120 ,
     label = "intervalo = 420-649 nm;
              PPFD = 337 um/m2/s;
              power= 74.4 W/m2",
     pos=3,
     cex = 0.9,  col="black")

#abline ( v=415, lty=1, col="navyblue", lwd=2 )
#abline ( v=584 , lty=1, col="navyblue", lwd=2 )


#text(x = 500, 
 #    y= 80 ,
 #    label = "BG.amb= 415-584 nm",
  #   pos=1,
   #  cex = 1.1,  col="navyblue")


#abline ( v=550, lty=1, col="darkred" , lwd=2)
#abline ( v=700 , lty=1, col="darkred", lwd=2 )

#text(x = 650, 
 #    y= 120 ,
  #   label = "Red.amb= 550-700 nm",
   #  pos=1,
    # cex = 1.1,  col="darkred")


###

plot (uW_nm.cm2 ~ lambda, 
      ylim = c(0,150),
      xlim=  c(400, 750),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c ( "power",unique(amb_GT.3$amb) ,sep="_",collapse = TRUE),
      data= amb_GT.3) 

abline ( v=420, lty=2, col="gray48" )
abline ( v=649 , lty=2, col="gray48" )
abline ( v=415, lty=1, col="navyblue", lwd=2 )
abline ( v=584 , lty=1, col="navyblue", lwd=2)

abline ( v=490, lty=1, col="black", lwd=2 )
abline ( v=560 , lty=1, col="black", lwd=2)


x1.pwr <- 606
x2.pwr <- 1086


#amb_GT.1$lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]

with (amb_GT.3, polygon(x=c(lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]),
                        y= c(0, uW_nm.cm2[x1.pwr:x2.pwr], 0),border ='navyblue',
                        col=scales::alpha( 'navyblue',.5)))

text(x = 510 , 
     y= 120 ,
     label = "inter = 415-584 nm;
              PPFD = 312.3 um/m2/s;
              power= 78.5 W/m2",
     #pos=2,
     cex = 0.9,  col="black")

####

plot (uW_nm.cm2 ~ lambda, 
      ylim = c(0,170),
      xlim=  c(400, 750),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c ( "power",unique(amb_GT.4$amb) ,sep="_",collapse = TRUE),
      data= amb_GT.4) 

abline ( v=420, lty=2, col="gray48" )
abline ( v=649 , lty=2, col="gray48" )
abline ( v=550, lty=1, col="darkred" , lwd=2)
abline ( v=700 , lty=1, col="darkred", lwd=2 )


abline ( v=490, lty=1, col="black", lwd=2 )
abline ( v=560 , lty=1, col="black", lwd=2)

x1.pwr <- 987
x2.pwr <- 1433


#amb_GT.1$lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]

with (amb_GT.4, polygon(x=c(lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]),
                        y= c(0, uW_nm.cm2[x1.pwr:x2.pwr], 0),border ='darkred',
                        col=scales::alpha( 'darkred',.5)))

text(x = 600 , 
     y= 120 ,
     label = "inter = 550-700 nm;
              PPFD = 396 um/m2/s;
              power= 74.7 W/m2",
     #pos=2,
     cex = 0.9,  col="black")

#####################################




text(x = 500, 
     y= 80 ,
     label = "BG.amb= 415-584 nm",
     pos=1,
     cex = 1.1,  col="navyblue")


abline ( v=550, lty=1, col="darkred" , lwd=2)
abline ( v=700 , lty=1, col="darkred", lwd=2 )

text(x = 650, 
     y= 120 ,
     label = "Red.amb= 550-700 nm",
     pos=1,
     cex = 1.1,  col="darkred")






x1=c(420:649, 649:420)
y1=c(0,rep (max(amb_GT.1$uW_nm.cm2),length(x1)-1))
length(x1)
polygon(x=x1, y=y1, col = "orange", 
        lty =1, lwd = 2, border = "red")


rect(x0=420,y0=max(amb_GT.1$uW_nm.cm2),x1=649,y1=y0, 
     col= rgb(0,0,1.0,alpha=0.5))

text(x = (dt.lmbd.x$l1 + dt.lmbd.x$l2)/2,
     y= max (dt.1.int.phot.total$uW_nm.cm2 - 10) ,
     label = str_c(dt.lmbd.x$l1 ,"_",dt.lmbd.x$l2),
     cex = 0.9,  col="gray48")


amb.luminico.GT.3




#encuentro indices de lambda = 400 y 700nm
rows_400nm <- which(grepl(400, dist_espectral_uw_raw$lambda))
rows_700nm <- which(grepl(700, dist_espectral_uw_raw$lambda))
index400nm <- rows_400nm[1]
index700nm <- rows_700nm[length(rows_700nm)]

#me quedo solo con el rango par
dist_espectral_uw <- dist_espectral_uw_raw[index400nm:index700nm,]

View (dist_espectral_uw)
str(dist_espectral_uw)


# Energía de un photon de una longitud de onda 
# E = hplanck*vluz/lambda
vluz    <- 3e8 #m/s
hplanck <- 6.626070e-34 #J.s
lambda_nm_a_m  <- 1e-9 # nm/m
uW_a_W  <- 1e-6 # uJ/J
Nav <- 6.02e23 # photones/uMol

E=
dist_espectral_variantes <- dist_espectral_uw %>%
                mutate (photones_nm.s.cm2 = uW_nm.cm2*uW_a_W
                							*lambda*lambda_nm_a_m
                							 /(hplanck* vluz)) %>%
                 mutate (umolphotones_nm.s.cm2 = photones_nm.s.cm2/Nav) 


str(dist_espectral_variantes)    

View(dist_espectral_variantes[index400nm: index700nm,])
 
# estructura de los datos
#svg (filename="./Figures/fig.X2.ambiente.luminico/dist.espectral.ps.sbf3.svg", 
 #    width=7, 
  #   height=5, 
   #  pointsize=12)
par(mfrow=c(2,1))
plot (uW_nm.cm2 ~ lambda,
      #ylim = c(0,250),
      xlim=  c(400,750),
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = "uW_nm.cm2",
      data= dist_espectral_variantes)

plot (photones_nm.s.cm2 ~ lambda,
      #ylim = c(0,250),
      xlim=  c(400,750),
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = "umolphotones_nm.s.cm2",
      data= dist_espectral_variantes)


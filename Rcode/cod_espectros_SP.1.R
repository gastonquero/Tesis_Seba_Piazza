##########################################################
#  Analisis de Datos Luz actinica PAM                    #
#   Ambientes luminicos JD.1                             #
#                                                        #
# Gaston Quero - Sebastian Fernandez - Sebastian Simondi #
#           10-11-2021                                   #
##########################################################

getwd ()
setwd ("R:/Tesis_Juan_Duhalde")

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

######### distribucion espectral  ######################

### seteo del directorio donde estan los datos

dir.data.base <- ("C:/Users/Usuario/OneDrive/Documentos/Paper_fotosintesis_dinamica")

#:\Users\Usuario\OneDrive\Documentos\Paper_fotosintesis_dinamica\Data\rawdata\Relevamiento_Luz\ambientes_luminicos_SP.1

dist_espectral_uw_SP.1_400_raw <- read.table (file =str_c(dir.data.base,"/Data/rawdata/Relevamiento_Luz/ambientes_luminicos_SP.1/sp.1.lb400_AbsoluteIrradiance_11-46-13-025.txt"),
                                              header = FALSE, sep = "\t",dec = ",", skip=30)

colnames (dist_espectral_uw_SP.1_400_raw) <- c("lambda", "uW_nm.cm2")


dist_espectral_uw_SP.1_400_raw <- dist_espectral_uw_SP.1_400_raw %>%
                                  dplyr::mutate (punto= "X")

#########  LB.800  #####################################


dist_espectral_uw_SP.1_800_raw <- read.table (file =str_c(dir.data.base,"/Data/rawdata/Relevamiento_Luz/ambientes_luminicos_SP.1/sp.1.lb800_AbsoluteIrradiance_11-49-39-524.txt"),
                                              
                                              header = FALSE, sep = "\t",dec = ",", skip=30)

colnames (dist_espectral_uw_SP.1_800_raw) <- c("lambda", "uW_nm.cm2")


dist_espectral_uw_SP.1_800_raw <- dist_espectral_uw_SP.1_800_raw %>%
                                  dplyr::mutate (punto= "X")


#########  LB.800a  #####################################

dist_espectral_uw_SP.1_800a_raw <- read.table (file =str_c(dir.data.base,"/Data/rawdata/Relevamiento_Luz/ambientes_luminicos_SP.1/SP.1_LB.800_AbsoluteIrradiance_1.txt"),
                                               header = FALSE, sep = "\t",dec = ",", skip=14)

colnames (dist_espectral_uw_SP.1_800a_raw) <- c("lambda", "uW_nm.cm2")

dist_espectral_uw_SP.1_800a_raw <- dist_espectral_uw_SP.1_800a_raw %>%
                                   dplyr::mutate (punto= "X")


#########  LB.Inv  #####################################

dist_espectral_uw_SP.1_Inv_raw <- read.table (file =str_c(dir.data.base,"/Data/rawdata/Relevamiento_Luz/ambientes_luminicos_SP.1/euca.invernaculo_AbsoluteIrradiance_13-00-51-880.txt"),
                                              header = FALSE, sep = "\t",dec = ",", skip=30)


colnames (dist_espectral_uw_SP.1_Inv_raw) <- c("lambda", "uW_nm.cm2")

dist_espectral_uw_SP.1_Inv_raw <- dist_espectral_uw_SP.1_Inv_raw %>%
                                  dplyr::mutate (punto= "X")








# Construccion de la matriz de lambdas
d1 <- tibble ( bandwidth = "D1", l1 = 400, l2= 425)
d2 <- tibble ( bandwidth = "D2", l1 = 425, l2= 490)
d3 <- tibble ( bandwidth = "D3", l1 = 490, l2= 560)
d4 <- tibble ( bandwidth = "D4", l1 = 560, l2= 585)
d5 <- tibble ( bandwidth = "D5", l1 = 585, l2= 640)
d6 <- tibble ( bandwidth = "D6", l1 = 640, l2= 700)
dPAR <- tibble ( bandwidth = "PAR", l1 = 400, l2= 700)
dintegr <- tibble ( bandwidth = "total", l1 = 380, l2= 780)

lambdas <- bind_rows (d1, d2, d3, d4, d5, d6,
                      dPAR,dintegr )


####
amb.luminico.LB.400 <- run_spectrum (dt = dist_espectral_uw_SP.1_400_raw, 
                                     lambdas = lambdas, id.survey="LB400", px="X" )

amb.luminico.LB.400  <- amb.luminico.LB.400 %>%
                        dplyr::mutate (relevamiento = "amb.luminico.SP.1.LB400")


amb.luminico.LB.800 <- run_spectrum (dt = dist_espectral_uw_SP.1_800_raw, 
                                     lambdas = lambdas, id.survey="LB800", px="X" )

amb.luminico.LB.800  <- amb.luminico.LB.800 %>%
                        dplyr::mutate (relevamiento = "amb.luminico.SP.1.LB800")

amb.luminico.LB.800a <- run_spectrum (dt = dist_espectral_uw_SP.1_800a_raw, 
                                     lambdas = lambdas, id.survey="LB800a", px="X" )

amb.luminico.LB.800a  <- amb.luminico.LB.800a %>%
                         dplyr::mutate (relevamiento = "amb.luminico.SP.1.LB800a")


amb.luminico.Inv <- run_spectrum (dt = dist_espectral_uw_SP.1_Inv_raw, 
                                  lambdas = lambdas, id.survey="Inv", px="X" )

amb.luminico.Inv  <- amb.luminico.Inv %>%
                     dplyr::mutate (relevamiento = "amb.luminico.SP.1.Inv")



calidades.luminicas.fot.SP1 <- bind_rows (amb.luminico.LB.400,
                                          amb.luminico.LB.800,
                                          amb.luminico.LB.800a,
                                          amb.luminico.Inv )


write_delim (calidades.luminicas.fot.SP1, file= str_c(dir.data.base,"/Data/procdata/calidades.luminicas.fot.SP1.txt"), delim = ",", na = "NA")

#### codigo para hacer la figura de los ambientes

amb_LB800  <- dist_espectral_uw_SP.1_800_raw  %>%
              dplyr::mutate (amb = "LB800")%>%
              dplyr::select (-punto)

amb_LB800a  <- dist_espectral_uw_SP.1_800a_raw  %>%
              dplyr::mutate (amb = "LB800a")%>%
              dplyr::select (-punto)

amb_LB400  <- dist_espectral_uw_SP.1_400_raw  %>%
              dplyr::mutate (amb = "LB400")%>%
              dplyr::select (-punto)

amb_Inv  <- dist_espectral_uw_SP.1_Inv_raw  %>%
            dplyr::mutate (amb = "INV")%>%
            dplyr::select (-punto)
             

#amb_GT <- bind_rows (amb_GT.1,amb_GT.2,amb_GT.3,amb_GT.4)             
             

plot (uW_nm.cm2 ~ lambda, 
      ylim = c(0,100),
      xlim=  c(400, 750),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c ( "power",unique(amb_LB800$amb) ,sep="_",collapse = TRUE),
      data= amb_LB800) 

abline ( v=400, lty=2, col="gray48" )
abline ( v=700 , lty=2, col="gray48" )


x1.pwr <- 565
x2.pwr <- 1433


amb_LB800$lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]

with (amb_LB800, polygon(x=c(lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]),
                       y= c(0, uW_nm.cm2[x1.pwr:x2.pwr], 0),border ='gray28',
                       col=scales::alpha( 'gray28',.3)))

text(x = 500 , 
     y= 80 ,
     label = "intervalo = 400-700 nm;
              PPFD = 750  um/m2/s;
              power= 16260 uW/cm2",
     pos=3,
     cex = 0.9,  col="black")

###

plot (uW_nm.cm2 ~ lambda, 
      ylim = c(0,130),
      xlim=  c(400, 750),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c ( "power",unique(amb_LB800a$amb) ,sep="_",collapse = TRUE),
      data= amb_LB800a) 

abline ( v=400, lty=2, col="gray48" )
abline ( v=700 , lty=2, col="gray48" )


x1.pwr <- 565
x2.pwr <- 1433


#amb_GT.1$lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]

with (amb_LB800a, polygon(x=c(lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]),
                        y= c(0, uW_nm.cm2[x1.pwr:x2.pwr], 0),border ='gray28',
                        col=scales::alpha( 'gray28',.3)))

text(x = 500 , 
     y= 100 ,
     label = "intervalo = 400-700 nm;
              PPFD = 818 um/m2/s;
              power= 16709.87 uW/cm2",
     pos=3,
     cex = 0.9,  col="black")


###

plot (uW_nm.cm2 ~ lambda, 
      ylim = c(0,120),
      xlim=  c(400, 750),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c ( "power",unique(amb_LB400$amb) ,sep="_",collapse = TRUE),
      data= amb_LB400) 


x1.pwr <- 565
x2.pwr <- 1433


#amb_GT.1$lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]

with (amb_LB400, polygon(x=c(lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]),
                        y= c(0, uW_nm.cm2[x1.pwr:x2.pwr], 0),border ='navyblue',
                        col=scales::alpha( 'navyblue',.5)))

text(x = 510 , 
     y= 100 ,
     label = "inter = 400-700 nm;
              PPFD = 413.3 um/m2/s;
              power= 8448.89 uW/cm2",
     #pos=2,
     cex = 0.9,  col="black")

####

plot (uW_nm.cm2 ~ lambda, 
      ylim = c(0,120),
      xlim=  c(400, 750),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c ( "power",unique(amb_Inv$amb) ,sep="_",collapse = TRUE),
      data= amb_Inv) 



abline ( v=400, lty=1, col="black", lwd=2 )
abline ( v=700 , lty=1, col="black", lwd=2)

x1.pwr <- 565
x2.pwr <- 1433


#amb_GT.1$lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]

with (amb_Inv, polygon(x=c(lambda [c(x1.pwr,x1.pwr:x2.pwr,x2.pwr)]),
                        y= c(0, uW_nm.cm2[x1.pwr:x2.pwr], 0),border ='darkred',
                        col=scales::alpha( 'darkred',.5)))

text(x = 500 , 
     y= 100 ,
     label = "inter = 400-700 nm;
              PPFD = 427.57 um/m2/s;
              power= 9389.69 uW/cm2",
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


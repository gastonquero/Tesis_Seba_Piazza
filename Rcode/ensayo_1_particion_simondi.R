##############################################################################
# Codigo para la funcion de PAM
# 
# Gaston Quero - Sebastian Simondi                                
# 20/01/2021
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

### Aca reingreso los datos la particion de simondi Fo


### E.1
E.1_part.simondi.fo <- read_delim (file ="./Data/procdata/phi.PAM.piazzas.Fo.txt", 
                            delim =",", quote = "\"",
                            escape_backslash = FALSE,
                            escape_double = TRUE, 
                            col_names = TRUE, 
                            col_types = NULL, 
                            locale = default_locale(), 
                            na = "NA")
unique(E.1_part.simondi.fo$sampling)

# genero las categoria de los ambientes 
E.1.inv <- E.1_part.simondi.fo %>%
           dplyr::filter ( sampling =="Inv" )


E.1.LB <- E.1_part.simondi.fo %>%
          dplyr::filter ( sampling !="Inv" ) %>%
          dplyr::mutate (x = pot) %>%
          tidyr::separate (x , c(NA,"X1",NA,NA, NA), sep ="_") %>%
          tidyr::separate (X1 , c("flujo", NA)) %>%
          dplyr::mutate (sampling = str_c (sampling, flujo, sep="_"))%>%
          dplyr::select (-flujo)

E.1.1 <-  bind_rows (E.1.inv , E.1.LB )

summary (E.1.1 )


#### reemplazo los valores negativos despues los vemos en detalle con Seba

E.1.1$phiPS2  [E.1.1$phiPS2  < 0]   <- NA

E.1.1$qp.fo  [E.1.1$qp.fo  < 0]   <- NA


E.1.1$phiNPQ [E.1.1$phiNPQ < 0]   <- NA


E.1.1$ppfd.act [E.1.1$ppfd.act  == 403]   <- 400

##### primera idea #####################
# modelar para tiempo final 
unique (E.1.1$pot)

t5.5  <- E.1.1 %>%
         dplyr::filter (time.min == 5.5) %>%
         dplyr::filter (pot != "LB_400.6_HHK1C3_3_Param.asc") %>%
         #dplyr::filter (ppfd.act < 1700 ) %>%
         dplyr::mutate (actinica = as.factor (ppfd.act ))%>%
         dplyr::mutate  (trat =(str_c (sampling , actinica , genotype )))

summary (t5.5)


t5.5$sampling <- as.factor (t5.5$sampling )


t5.5$genotype <- as.factor (t5.5$genotype)



# modelo completo del invernaculo 
## actinica categorica, 
t5.5.inv <- t5.5  %>%
            dplyr::filter ( sampling == "Inv")


## model 
lm.phiPS2.inv <- lm (phiPS2 ~  actinica + genotype + 
                      actinica * genotype 
                    , data= t5.5.inv )

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiPS2.inv$fitted.values, lm.phiPS2.inv$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiPS2.inv")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiPS2.inv$residuals, main="phiPS2.inv",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiPS2.inv$residuals),
              sd (lm.phiPS2.inv$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiPS2.inv$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiPS2.inv")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiPS2.inv$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiPS2 ~ trat, data = t5.5.inv)

# ANAVA
(anava.lm.phiPS2.inv  <- anova (lm.phiPS2.inv))
my_comparisons <- list( c("200", "400"), c("850", "400"), c("1700", "400"), c("850", "1700"))
p.inv <- ggboxplot (t5.5.inv, x = "actinica", y = "phiPS2",
                ylim=c(0,1),
                title = "Inv.",
               color = "genotype",
               palette =c("navyblue", "orange"),
               add = c("mean_sd","jitter"),
               add.params = list(color = "genotype")) + 
               stat_compare_means (comparisons = my_comparisons, 
                        label = "p.signif", method = "t.test")

ggexport (p.inv, filename = "./Figures/phiPS2.inv.tiff",
          width = 700, height = 500)



# modelo completo del LB_400
## actinica categorica,

t5.5.LB400 <- t5.5  %>%
              dplyr::filter ( sampling == "LB_400")


## model 
lm.phiPS2.LB400 <- lm (phiPS2 ~  actinica + genotype + 
                       actinica * genotype 
                     , data= t5.5.LB400 )

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiPS2.LB400$fitted.values, lm.phiPS2.LB400$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiPS2.LB400")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiPS2.LB400$residuals, main="phiPS2.LB400",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiPS2.LB400$residuals),
              sd (lm.phiPS2.LB400$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiPS2.LB400$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiPS2.LB400")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiPS2.LB400$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiPS2 ~ trat, data = t5.5.LB400)

# ANAVA
(anava.lm.phiPS2.LB400  <- anova (lm.phiPS2.LB400))

my_comparisons <- list( c("200", "400"), c("850", "400"), c("1700", "400"), c("850", "1700"))

p.400 <- ggboxplot (t5.5.LB400, x = "actinica", y = "phiPS2",
                    title = "LB_400",
                ylim=c(0,1),
                color = "genotype",
                palette =c("navyblue", "orange"),
                add = c("mean_sd","jitter"),
                add.params = list(color = "genotype")) + 
           stat_compare_means (comparisons = my_comparisons, 
                        label = "p.signif", method = "t.test")

ggexport (p.400, filename = "./Figures/phiPS2.LB400.tiff",
          width = 700, height = 500)

# modelo completo del LB_800
## actinica categorica,

t5.5.LB800 <- t5.5  %>%
              dplyr::filter ( sampling == "LB_800") %>%
              dplyr::filter (actinica != "1700")

summary (t5.5.LB800 )

## model 800
lm.phiPS2.LB800 <- lm (phiPS2 ~  actinica + genotype + 
                         actinica * genotype 
                       , data= t5.5.LB800 )

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiPS2.LB800$fitted.values, lm.phiPS2.LB800$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiPS2.LB800")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiPS2.LB800$residuals, main="phiPS2.LB800",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiPS2.LB800$residuals),
              sd (lm.phiPS2.LB800$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiPS2.LB800$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiPS2.LB800")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiPS2.LB800$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiPS2 ~ trat, data = t5.5.LB800)

# ANAVA
(anava.lm.phiPS2.LB800  <- anova (lm.phiPS2.LB800))


my_comparisons <- list( c("200", "400"), c("850", "400"), c("1700", "400"), c("850", "1700"))

p.800 <- ggboxplot (t5.5.LB800, x = "actinica", y = "phiPS2",
                    title = "LB_800",
                    ylim=c(0,1),
                    color = "genotype",
                    palette =c("navyblue", "orange"),
                    add = c("mean_sd","jitter"),
                    add.params = list(color = "genotype")) + 
                    stat_compare_means (comparisons = my_comparisons, 
                            label = "p.signif", method = "t.test")



my_comparisons.800 <- list( c("200", "400"), c("850", "400"))

p.800.1 <- ggboxplot (t5.5.LB800, x = "actinica", y = "phiPS2",
                    title = "LB_800",
                    ylim=c(0,1),
                    color = "genotype",
                    palette =c("navyblue", "orange"),
                    add = c("mean_sd","jitter"),
                    add.params = list(color = "genotype")) + 
  stat_compare_means (comparisons = my_comparisons.800, 
                      label = "p.signif",method = "t.test")



ggexport (p.800, filename = "./Figures/phiPS2.LB800.tiff",
          width = 700, height = 500)

ggexport (p.800.1, filename = "./Figures/phiPS2.LB800.1.tiff",
          width = 700, height = 500)

# modelo conjunto 
# modelo completo del LB_800
## actinica categorica,

t5.5.LB <- t5.5  %>%
  dplyr::filter ( sampling != "Inv") #%>%
  #dplyr::filter (actinica != "1700")

summary (t5.5.LB )

## model LB
lm.phiPS2.LB <- lm (phiPS2 ~  sampling + actinica + genotype + 
                              sampling * actinica +
                              sampling * genotype +
                              actinica * genotype +
                              sampling * actinica * genotype
                       , data= t5.5.LB )

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiPS2.LB$fitted.values, lm.phiPS2.LB$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiPS2.LB")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiPS2.LB$residuals, main="phiPS2.LB",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiPS2.LB$residuals),
              sd (lm.phiPS2.LB$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiPS2.LB$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiPS2.LB")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiPS2.LB$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiPS2 ~ trat, data = t5.5.LB)

# ANAVA
(anava.lm.phiPS2.LB  <- anova (lm.phiPS2.LB))

t5.5.LB.1C3 <- t5.5.LB %>%
               dplyr::filter (genotype == "HHK1C3")



p.1C3 <- ggboxplot (t5.5.LB.1C3, x = "actinica", y = "phiPS2", facet.by = "sampling",
        ylim=c(0,1),
        title = "HHK1C3",
        color = "sampling", 
        palette =c( "red", "green"),
        add = c("mean_sd","jitter"),
        add.params = list(color = "sampling")) + 
  stat_compare_means (comparisons = my_comparisons, 
                      label = "p.signif",method = "t.test")

t5.5.LB.2C5 <- t5.5.LB %>%
  dplyr::filter (genotype == "HHK2C5")



p.2C5 <-ggboxplot (t5.5.LB.2C5, x = "actinica", y = "phiPS2", facet.by = "sampling",
           ylim=c(0,1),
           title = "HHK2C5",
           color = "sampling", 
           palette =c( "red", "green"),
           add = c("mean_sd","jitter"),
           add.params = list(color = "sampling")) + 
  stat_compare_means (comparisons = my_comparisons, 
                      label = "p.signif",method = "t.test")


ggarrange (p.1C3, p.2C5, 
           ncol=1,
           nrow = 2)



ggline (t5.5.LB, x = "trat", y = "len", add = "mean_se",
       color = "supp", palette = "jco")+
  stat_compare_means(aes(group = supp), label = "p.signif", 
                     label.y = c(16, 25, 29))






# Visualize the expression profile
ggboxplot  (t5.5.LB, x = "trat",  y = "phiPS2", color = "trat", 
            add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)


+
  geom_hline(yintercept = mean(myeloma$DEPDC1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 1600)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all





my_comparisons <- list( c("200", "400"), c("850", "400"), c("1700", "400"), c("850", "1700"))

p.800 <- ggboxplot (t5.5.LB, x = "actinica", y = "phiPS2",
                    title = "LB_800",
                    ylim=c(0,1),
                    color = "genotype",
                    palette =c("navyblue", "orange"),
                    add = c("mean_sd","jitter"),
                    add.params = list(color = "genotype")) + 
  stat_compare_means (comparisons = my_comparisons, 
                      label = "p.signif", method = "t.test")



my_comparisons.800 <- list( c("200", "400"), c("850", "400"))

p.800.1 <- ggboxplot (t5.5.LB, x = "actinica", y = "phiPS2",
                      title = "LB_800",
                      ylim=c(0,1),
                      color = "genotype",
                      palette =c("navyblue", "orange"),
                      add = c("mean_sd","jitter"),
                      add.params = list(color = "genotype")) + 
  stat_compare_means (comparisons = my_comparisons.800, 
                      label = "p.signif",method = "t.test")



ggexport (p.800, filename = "./Figures/phiPS2.LB.tiff",
          width = 700, height = 500)

ggexport (p.800.1, filename = "./Figures/phiPS2.LB.1.tiff",
          width = 700, height = 500)





p + stat_compare_means(aes(group = genotype),label = "p.signif")

p + stat_compare_means (comparisons = my_comparisons, 
                       label = "p.signif", method = "t.test")

p + stat_compare_means(method = "anova", label.y = 0.7) +      # Add global p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "400")           




                        
                        
                        
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value





### analizamos el efecto principal de la luz actinica 

### 200
em.phiPS2.200.inv <- emmeans (lm.phiPS2.inv, ~ genotype, 
                                at = list (actinica = "200"))

(cr.phiPS2.200.inv <- contrast (em.phiPS2.200.inv, method = "pairwise"))

cld.phiPS2.HHK1C3.LB <- cld (em.phiPS2.HHK1C3.LB, sort=FALSE)

blue.cld.phiPS2.HHK1C3.LB <- cbind(genotipo= "HHK1C3", parametro ="phiPS2",cld.phiPS2.HHK1C3.LB )




unique (t5.5$pot )

t5.5.inv <-  t5.5 %>%
             dplyr::filter (sampling == "Inv")


ggscatter (t5.5.inv  , x = "ppfd.act", y = "phiPS2",
           add = "loess", size = 3,
           facet.by = "genotype" ,
           color = "pot",
           palette = "jco")


ggscatter (t5.5.inv  , x = "ppfd.act", y = "phiPS2",
           add = "loess", size = 3,
           facet.by = "genotype" ,
           color = "genotype",
           palette = "jco")




t5.5.LB400 <- t5.5 %>%
              dplyr::filter (time.min == 5.5) %>%
              dplyr::filter (sampling == "LB_400")


ggscatter (t5.5.LB400  , x = "ppfd.act", y = "phiPS2",
           add = "loess", size = 3,
           facet.by = "genotype" ,
           color = "pot",
           palette = "jco")


ggscatter (t5.5.LB400  , x = "ppfd.act", y = "phiPS2",
           add = "loess", size = 3,
           facet.by = "genotype" ,
           color = "genotype",
           palette = "jco")




t5.5.inv.400 <- t5.5 %>%
                dplyr::filter (sampling != "LB_800")

ggscatter (t5.5.inv.400 , x = "ppfd.act", y = "phiPS2",
           add = "loess", size = 3,
           facet.by = "sampling" ,
           color = "genotype",
           palette = "jco")


ggscatter (t5.5.inv.400 , x = "ppfd.act", y = "phiPS2",
           add = "loess", size = 3,
           facet.by = "genotype",
           color = "sampling",
           palette = "jco")



ggscatter (t5.5.LB400  , x = "ppfd.act", y = "phiPS2",
           add = "loess", size = 3,
           facet.by = "genotype",
           color = "sampling",
           palette = c("black", "red", "green"))


### por pot 

unique (t5.5$pot)

#t5.5 <- t5.5 %>%
 #       dplyr::filter (pot!= "LB_800.5_HHK2C5_3_Param.asc") 

pot.x <- t5.5 %>%
         dplyr::filter (sampling == "LB_400" ) 

ggscatter (pot.x  , x = "ppfd.act", y = "phiPS2", 
           facet.by = "genotype",
           ylim =c(0.15, 0.8),
           add = "loess", size = 3,
           color="pot",
           palette ="jco",
           )

pot.y <- t5.5 %>%
  dplyr::filter (sampling == "LB_800" ) 


ggscatter (pot.y  , x = "ppfd.act", y = "phiPS2", 
           facet.by = "genotype",
           ylim =c(0.15, 0.8),
           add = "loess", size = 3,
           color="pot",
           palette ="jco",
)



#### solo LB
t5.5.LB <- t5.5 %>%
           dplyr::filter (sampling != "Inv" )

t5.5.LB$sampling <- as.factor (t5.5.LB$sampling )

t5.5.LB$ppfd.act <- as.factor (t5.5.LB$ppfd.act )

t5.5.LB$genotype <- as.factor (t5.5.LB$genotype)


t5.5.LB <- t5.5.LB %>%
            dplyr::mutate (trat = str_c (sampling ,ppfd.act,genotype  ))
 

str (t5.5.LB )
summary (t5.5.LB)

#### vamos a modelar phiPS2 ####

## actinica categorica, 

## model 
lm.phiPS2.LB <- lm (phiPS2 ~ sampling + ppfd.act + genotype + 
                   sampling * ppfd.act + 
                   sampling * genotype + 
                   ppfd.act * genotype +
                   sampling * ppfd.act * genotype
                   , data= t5.5.LB )

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiPS2.LB$fitted.values, lm.phiPS2.LB$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiPS2.LB")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiPS2.LB$residuals, main="phiPS2.LB",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiPS2.LB$residuals),
              sd (lm.phiPS2.LB$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiPS2.LB$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiPS2.LB")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiPS2.LB$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiPS2 ~ trat, data = t5.5.LB)

# ANAVA
(anava.lm.phiPS2  <- anova (lm.phiPS2.LB))

### analizamos al interaccion sampling_genotype

unique (t5.5$genotype)

### HHK1C3
em.phiPS2.HHK1C3.LB <- emmeans (lm.phiPS2.LB, ~ sampling:ppfd.act , 
                               at = list (genotype = "HHK1C3"))

(cr.phiPS2.HHK1C3.LB <- contrast (em.phiPS2.HHK1C3.LB , method = "pairwise"))

cld.phiPS2.HHK1C3.LB <- cld (em.phiPS2.HHK1C3.LB, sort=FALSE)

blue.cld.phiPS2.HHK1C3.LB <- cbind(genotipo= "HHK1C3", parametro ="phiPS2",cld.phiPS2.HHK1C3.LB )

#### 400 
phiPS2.ref.1C3.400 <- blue.cld.phiPS2.HHK1C3.LB %>%
                         dplyr::filter (sampling == "LB_400" ) %>%
                         dplyr::filter (ppfd.act == "403" ) %>% 
                         dplyr::select (emmean )

phiPS2.HHK1C3.LB.400 <- blue.cld.phiPS2.HHK1C3.LB %>%
                 dplyr::filter (sampling == "LB_400" ) %>%
                 dplyr::select (genotipo,parametro, sampling, ppfd.act,emmean ) %>%
                 dplyr::mutate (phi.relative = emmean/phiPS2.ref.1C3.400$emmean)

phiPS2.HHK1C3.LB.400$ppfd.act <- factor (phiPS2.HHK1C3.LB.400$ppfd.act, levels = c("200", "403", "850", "1700" ))


### p1C3.400.re
p1C3.400.rel <- ggdotplot (phiPS2.HHK1C3.LB.400, x ="ppfd.act", y = "phi.relative",
           ylim=c(0,2.25),
           size=2.5,
           title = "p1C3.400",
           ylab = "phiPSII.relativo",
           color = "black", fill = "black") + 
            geom_hline (yintercept = 1, lty=2) + 
            geom_vline (xintercept = "403", lty=2)

### p1C3.400.abs
p1C3.400.abs <- ggdotplot (phiPS2.HHK1C3.LB.400, x ="ppfd.act", y = "emmean",
                         ylim=c(0,0.8),
                         size=2,
                         title = "p1C3.400",
                         ylab = "phiPSII.abs",
                         color = "black", fill = "black") + 
  geom_hline (yintercept =phiPS2.ref.1C3.400$emmean, lty=2) + 
  geom_vline (xintercept = "403", lty=2)



### 800

phiPS2.ref.1C3.800 <- blue.cld.phiPS2.HHK1C3.LB %>%
                      dplyr::filter (sampling == "LB_800" ) %>%
                      dplyr::filter (ppfd.act == "850" ) %>% 
                       dplyr::select (emmean )

phiPS2.HHK1C3.LB.800 <- blue.cld.phiPS2.HHK1C3.LB %>%
  dplyr::filter (sampling == "LB_800" ) %>%
  dplyr::select (genotipo,parametro, sampling, ppfd.act,emmean ) %>%
  dplyr::mutate (phi.relative = emmean/phiPS2.ref.1C3.800$emmean)

phiPS2.HHK1C3.LB.800$ppfd.act <- factor (phiPS2.HHK1C3.LB.800$ppfd.act, levels = c("200", "403", "850", "1700" ))


# plot p1C3.800.re
p1C3.800.rel <- ggdotplot (phiPS2.HHK1C3.LB.800, x ="ppfd.act", y = "phi.relative",
                       ylim=c(0,2.25),
                       size=2,
                       title = "p1C3.800",
                       ylab = "phiPSII.relativo",
                       color = "black", fill = "black") + 
  geom_hline (yintercept = 1, lty=2) + 
  geom_vline (xintercept = "850", lty=2)

# plot p1C3.800.abs
p1C3.800.abs <- ggdotplot (phiPS2.HHK1C3.LB.800, x ="ppfd.act", y = "emmean",
                           ylim=c(0,0.8),
                           size=2,
                         title = "p1C3.800",
                         ylab = "phiPSII.abs",
                         color = "black", fill = "black") + 
  geom_hline (yintercept =phiPS2.ref.1C3.800$emmean, lty=2) + 
  geom_vline (xintercept = "850", lty=2)


### 2C5 
### HHK2C5
em.phiPS2.HHK2C5.LB <- emmeans (lm.phiPS2.LB, ~ sampling:ppfd.act , 
                                at = list (genotype = "HHK2C5"))

(cr.phiPS2.HHK2C5.LB <- contrast (em.phiPS2.HHK2C5.LB , method = "pairwise"))

cld.phiPS2.HHK2C5.LB <- cld (em.phiPS2.HHK2C5.LB, sort=FALSE)

blue.cld.phiPS2.HHK2C5.LB <- cbind(genotipo= "HHK2C5", parametro ="phiPS2",cld.phiPS2.HHK2C5.LB )

#### 400 
phiPS2.ref.2C5.400 <- blue.cld.phiPS2.HHK2C5.LB %>%
  dplyr::filter (sampling == "LB_400" ) %>%
  dplyr::filter (ppfd.act == "403" ) %>% 
  dplyr::select (emmean )

phiPS2.HHK2C5.LB.400 <- blue.cld.phiPS2.HHK2C5.LB %>%
  dplyr::filter (sampling == "LB_400" ) %>%
  dplyr::select (genotipo,parametro, sampling, ppfd.act,emmean ) %>%
  dplyr::mutate (phi.relative = emmean/phiPS2.ref.2C5.400$emmean)

phiPS2.HHK2C5.LB.400$ppfd.act <- factor (phiPS2.HHK2C5.LB.400$ppfd.act, levels = c("200", "403", "850", "1700" ))


# plot p2C5.400.rel
p2C5.400.rel <- ggdotplot (phiPS2.HHK2C5.LB.400, x ="ppfd.act", y = "phi.relative",
                           ylim=c(0,2.25),
                           size=2,
                       title = "p2C5.400",
                       ylab = "phiPSII.relativo",
                       color = "black", fill = "black") + 
  geom_hline (yintercept = 1, lty=2) + 
  geom_vline (xintercept = "403", lty=2)

# plot p2C5.400.abs
p2C5.400.abs <- ggdotplot (phiPS2.HHK2C5.LB.400, x ="ppfd.act", y = "emmean",
                           ylim=c(0,0.8),
                           size=2,
                           title = "p2C5.400",
                           ylab = "phiPSII.abs",
                           color = "black", fill = "black") + 
  geom_hline (yintercept =phiPS2.ref.2C5.400$emmean, lty=2) + 
  geom_vline (xintercept = "403", lty=2)




### 800

phiPS2.ref.2C5.800 <- blue.cld.phiPS2.HHK2C5.LB %>%
  dplyr::filter (sampling == "LB_800" ) %>%
  dplyr::filter (ppfd.act == "850" ) %>% 
  dplyr::select (emmean )

phiPS2.HHK2C5.LB.800 <- blue.cld.phiPS2.HHK2C5.LB %>%
  dplyr::filter (sampling == "LB_800" ) %>%
  dplyr::select (genotipo,parametro, sampling, ppfd.act,emmean ) %>%
  dplyr::mutate (phi.relative = emmean/phiPS2.ref.2C5.800$emmean)

phiPS2.HHK2C5.LB.800$ppfd.act <- factor (phiPS2.HHK2C5.LB.800$ppfd.act, levels = c("200", "403", "850", "1700" ))

#plot p2C5.800.rel
p2C5.800.rel <- ggdotplot (phiPS2.HHK2C5.LB.800, x ="ppfd.act", y = "phi.relative",
                           ylim=c(0,2.25),
                           size=2,
                       title = "p2C5.800",
                       ylab = "phiPSII.relativo",
                       color = "black", fill = "black") + 
  geom_hline (yintercept = 1, lty=2) + 
  geom_vline (xintercept = "850", lty=2)


#plot p2C5.800.abs

p2C5.800.abs <- ggdotplot (phiPS2.HHK2C5.LB.800, x ="ppfd.act", y = "emmean",
                           ylim=c(0,0.8),
                           size=2,
                           title = "p2C5.800",
                           ylab = "phiPSII.abs",
                           color = "black", fill = "black") + 
  geom_hline (yintercept =phiPS2.ref.2C5.800$emmean, lty=2) + 
  geom_vline (xintercept = "850", lty=2)



blue.phiPS2 <- bind_rows (blue.cld.phiPS2.HHK1C3.LB, blue.cld.phiPS2.HHK2C5.LB )

write_delim(
  blue.phiPS2,
  file= "./Data/procdata/blue.phiPS2.txt",
  delim = ";",
  na = "NA",
  append = FALSE)

write_excel_csv2(
  blue.phiPS2,
  file= "./Data/procdata/blue.phiPS2.csv"
)

ggarrange (p1C3.400.rel, p1C3.800.rel,
           p2C5.400.rel, p2C5.800.rel, 
           ncol = 2,
           nrow = 2)


ggarrange (p1C3.400.abs, p1C3.800.abs,
           p2C5.400.abs, p2C5.800.abs, 
           ncol = 2,
           nrow = 2)

##### actinica numerica

blue.phiPS2.200 <- blue.phiPS2 %>%
                   dplyr::filter (ppfd.act == "200" ) %>%
                   dplyr::mutate (ppfd.act.num = 200)


blue.phiPS2.400 <- blue.phiPS2 %>%
                   dplyr::filter (ppfd.act == "403" ) %>%
                   dplyr::mutate (ppfd.act.num = 400)


blue.phiPS2.850 <- blue.phiPS2 %>%
                   dplyr::filter (ppfd.act == "850" ) %>%
                   dplyr::mutate (ppfd.act.num = 850)


blue.phiPS2.1700 <- blue.phiPS2 %>%
                    dplyr::filter (ppfd.act == "1700" ) %>%
                    dplyr::mutate (ppfd.act.num = 1700)


blue.phiPS2.num <- bind_rows (blue.phiPS2.200,
                              blue.phiPS2.400,
                              blue.phiPS2.850,
                              blue.phiPS2.1700)



ggscatter (blue.phiPS2.num, x = "ppfd.act.num", y = "emmean",
           add = "loess", size = 3,
           facet.by = "sampling",
          color = "genotipo",
          palette = c("navyblue", "orange"))





p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(aes(color = gear))

# Default plot
p

# Use theme_pubr()
p + theme_pubr()

# Format labels
p + labs_pubr()

### HHK2C5
em.phiPS2.HHK2C5.LB <- emmeans (lm.phiPS2.LB, ~ sampling:ppfd.act , 
                                at = list (genotype = "HHK2C5"))

(cr.phiPS2.HHK2C5.LB <- contrast (em.phiPS2.HHK2C5.LB , method = "pairwise"))

cld.phiPS2.HHK2C5.LB <- cld (em.phiPS2.HHK2C5.LB, sort=FALSE)

blue.cld.phiPS2.HHK2C5.LB <- cbind(genotipo= "HHK2C5", parametro ="phiPS2",cld.phiPS2.HHK2C5.LB )







phiPS2.geno <- ggline (t5.5, x = "ppfd.act", y = "phiPS2", 
                       ylim=c(0,1),point.size=2,
                       facet.by = "genotype",
                       desc_stat = "mean_sd",
                       fill = "sampling", color = "sampling",
                       add = c("mean_sd"), palette = c("gray48", "orange", "darkred"))

ggexport (phiPS2.geno, filename = "./Figures/E1.phiPS2.geno.tiff",
          width = 700, height = 500)



ggline (t5.5, x = "ppfd.act", y = "phiPS2", 
        ylim=c(0,1),point.size=2,
        facet.by = "genotype",
        desc_stat = "mean_sd",
        fill = "sampling", color = "sampling",
        add = c("mean_sd"), palette = c("gray48", "orange", "darkred"))




qp.fo.geno <- ggline (t5.5, x = "ppfd.act", y = "qp.fo", 
                       ylim=c(0,1),point.size=2,
                       facet.by = "genotype",
                       desc_stat = "mean_sd",
                       fill = "sampling", color = "sampling",
                       add = c("mean_sd"), palette = c("gray48", "orange", "darkred"))

ggexport (qp.fo.geno , filename = "./Figures/E1.qp.fo.geno.tiff",
          width = 700, height = 500)




sp.fo.geno  <- ggline (t5.5, x = "ppfd.act", y = "sp.fo", 
                 ylim=c(0,1),point.size=2,
                 facet.by = "genotype",
                 desc_stat = "mean_sd",
                 fill = "sampling", color = "sampling",
                 add = c("mean_sd"), palette = c("gray48", "orange", "darkred"))


ggexport (sp.fo.geno , filename = "./Figures/E1.sp.fo.geno.tiff",
          width = 700, height = 500)

unique (t5.5$genotype)
#### relacion 
t5.5.2C5 <- t5.5 %>%
            dplyr::filter ( genotype == "HHK2C5")

ggscatter (t5.5.2C5 , x = "qp.fo" , y = "phiPS2", 
          ylim=c(0,1),point.size=2,
          facet.by = "ppfd.act",
          fill = "sampling", color = "sampling",
          #add = "reg.line",
          palette = c("gray48", "orange", "darkred"))
          
          
ggscatter (t5.5.2C5 , x = "sp.fo" , y = "phiPS2", 
           ylim=c(0,1),point.size=2,
           facet.by = "ppfd.act",
           fill = "sampling", color = "sampling",
           #add = "reg.line",
           palette = c("gray48", "orange", "darkred"))


######
t5.5.1C3 <- t5.5 %>%
            dplyr::filter ( genotype == "HHK1C3")

ggscatter (t5.5.1C3 , x = "qp.fo" , y = "phiPS2", 
           ylim=c(0,1),point.size=2,
           facet.by = "ppfd.act",
           fill = "sampling", color = "sampling",
           #add = "reg.line",
           palette = c("gray48", "orange", "darkred"))


ggscatter (t5.5.1C3 , x = "sp.fo" , y = "phiPS2", 
           ylim=c(0,1),point.size=2,
           facet.by = "ppfd.act",
           fill = "sampling", color = "sampling",
           #add = "reg.line",
           palette = c("gray48", "orange", "darkred"))


phiNPQ.geno <- ggline (t5.5, x = "ppfd.act", y = "phiNPQ", 
                       ylim=c(0,1),point.size=2,
                       facet.by = "genotype",
                       desc_stat = "mean_sd",
                       fill = "sampling", color = "sampling",
                       add = c("mean_sd"), palette = c("gray48", "orange", "darkred"))

ggexport (phiNPQ.geno, filename = "./Figures/E1.phiNPQ.geno.tiff",
          width = 700, height = 500)





phiNO.geno <- ggline (t5.5, x = "ppfd.act", y = "phiNO", 
                      ylim=c(0,1),point.size=2,
                      facet.by = "genotype",
                      desc_stat = "mean_sd",
                      fill = "sampling", color = "sampling",
                      add = c("mean_sd"), palette = c("gray48", "orange", "darkred"))

ggexport (phiNO.geno, filename = "./Figures/E1.phiNO.geno.tiff",
          width = 700, height = 500)




phiPS2.sampling  <- ggline (t5.5, x = "ppfd.act", y = "phiPS2", 
                            ylim=c(0,1),
                            point.size=2,
                            facet.by = "sampling",
                            desc_stat = "mean_sd",
                            fill = "genotype", color = "genotype",
                            add = c("mean_sd"), palette = c("navyblue", "darkorange"))

ggexport (phiPS2.sampling, filename = "./Figures/E1.phiPS2.sampling.tiff",
          width = 700, height = 500)


phiNPQ.sampling  <- ggline (t5.5, x = "ppfd.act", y = "phiNPQ", 
                            ylim=c(0,1),
                            point.size=2,
                            facet.by = "sampling",
                            desc_stat = "mean_sd",
                            fill = "genotype", color = "genotype",
                            add = c("mean_sd"), palette = c("navyblue", "darkorange"))

ggexport (phiNPQ.sampling, filename = "./Figures/E1.phiNPQ.sampling.tiff",
          width = 700, height = 500)


phiNO.sampling  <- ggline (t5.5, x = "ppfd.act", y = "phiNO", 
                           ylim=c(0,1),
                           point.size=2,
                           facet.by = "sampling",
                           desc_stat = "mean_sd",
                           fill = "genotype", color = "genotype",
                           add = c("mean_sd"), palette = c("navyblue", "darkorange"))

ggexport (phiNO.sampling, filename = "./Figures/E1.phiNO.sampling.tiff",
          width = 700, height = 500)


# idea del tiempo #########

E.1.1$ppfd.act <- as.character (E.1.1$ppfd.act)




unique (E.1.1$sampling )

### 
Inv <- E.1.1 %>%
       dplyr::filter (sampling == "Inv") 

Inv.H1C3 <- Inv  %>%
            dplyr::filter ( genotype == "HHK1C3")

unique (Inv.H1C3$ppfd.act )

Inv.H1C3$ppfd.act <- factor (Inv.H1C3$ppfd.act, levels = c("200", "403", "850","1700" ))

plot.1.Inv.H1C3 <- ggscatter (Inv.H1C3 , x = "sp.fo" , y = "phiPS2", 
                  title = "Inv.H1C3",
                 ylim=c(0,1),point.size=2,
                 cor.coef =TRUE,
                 cor.coeff.args = list(method = "pearson"),
                 facet.by = "ppfd.act",
                 fill = "ppfd.act", color = "ppfd.act",
                 add = "reg.line",
                 palette = c("orange", "red", "blue", "navyblue"))

ggexport (plot.1.Inv.H1C3 , filename = "./Figures/E1.plot.1.Inv.H1C3.tiff",
          width = 700, height = 500)


plot.2.Inv.H1C3 <- ggscatter (Inv.H1C3 , x = "qp.fo" , y = "phiPS2", 
                              title = "Inv.H1C3",
                              ylim=c(0,1),point.size=2,
                              cor.coef =TRUE,
                              cor.coeff.args = list(method = "pearson"),
                              facet.by = "ppfd.act",
                              fill = "ppfd.act", color = "ppfd.act",
                              add = "reg.line",
                              palette = c("orange", "red", "blue", "navyblue"))

ggexport (plot.2.Inv.H1C3 , filename = "./Figures/E1.plot.2.Inv.H1C3.tiff",
          width = 700, height = 500)



Inv.H2C5 <- Inv  %>%
            dplyr::filter ( genotype == "HHK2C5")

unique (Inv.H2C5$ppfd.act )

Inv.H2C5$ppfd.act <- factor (Inv.H2C5$ppfd.act, levels = c("200", "403", "850","1700" ))

plot.1.Inv.H2C5 <- ggscatter (Inv.H2C5 , x = "sp.fo" , y = "phiPS2", 
                   title = "Inv.H2C5",
                   ylim=c(0,1),point.size=2,
                   cor.coef =TRUE,
                   cor.coeff.args = list(method = "pearson"),
                   facet.by = "ppfd.act",
                   fill = "ppfd.act", color = "ppfd.act",
                   add = "reg.line",
                   palette = c("orange", "red", "blue", "navyblue"))


ggexport (plot.1.Inv.H2C5 , filename = "./Figures/E1.plot.1.Inv.H2C5.tiff",
          width = 700, height = 500)



plot.2.Inv.H2C5 <- ggscatter (Inv.H2C5 , x = "qp.fo" , y = "phiPS2", 
                  title = "Inv.H2C5",
                  ylim=c(0,1),point.size=2,
                  cor.coef =TRUE,
                  cor.coeff.args = list(method = "pearson"),
                  facet.by = "ppfd.act",
                  fill = "ppfd.act", color = "ppfd.act",
                  add = "reg.line",
                  palette = c("orange", "red", "blue", "navyblue"))

ggexport (plot.2.Inv.H2C5 , filename = "./Figures/E1.plot.2.Inv.H2C5.tiff",
          width = 700, height = 500)


### 
LB_400 <- E.1.1 %>%
          dplyr::filter (sampling == "LB_400") 

LB_400.H1C3 <- LB_400  %>%
               dplyr::filter ( genotype == "HHK1C3")

unique (LB_400.H1C3$ppfd.act )

LB_400.H1C3$ppfd.act <- factor (LB_400.H1C3$ppfd.act, levels = c("200", "403", "850","1700" ))

plot.1.LB_400.H1C3 <- ggscatter (LB_400.H1C3 , x = "sp.fo" , y = "phiPS2", 
                      title = "LB_400.H1C3",
                      ylim=c(0,1),point.size=2,
                      cor.coef =TRUE,
                      cor.coeff.args = list(method = "pearson"),
                      facet.by = "ppfd.act",
                      fill = "ppfd.act", color = "ppfd.act",
                      add = "reg.line",
                      palette = c("orange", "red", "blue", "navyblue"))


ggexport (plot.1.LB_400.H1C3, filename = "./Figures/E1.plot.1.LB_400.H1C3.tiff",
          width = 700, height = 500)



plot.2.LB_400.H1C3 <- ggscatter (LB_400.H1C3 , x = "qp.fo" , y = "phiPS2", 
                      title = "LB_400.H1C3",
                      ylim=c(0,1),
                      point.size=2,
                      cor.coef =TRUE,
                      cor.coeff.args = list(method = "pearson"),
                      facet.by = "ppfd.act",
                      fill = "ppfd.act", color = "ppfd.act",
                      add = "reg.line",
                      palette = c("orange", "red", "blue", "navyblue"))


ggexport (plot.2.LB_400.H1C3, filename = "./Figures/E1.plot.2.LB_400.H1C3.tiff",
          width = 700, height = 500)





LB_400.H2C5 <- LB_400  %>%
  dplyr::filter ( genotype == "HHK2C5")

unique (LB_400.H2C5$ppfd.act )

LB_400.H2C5$ppfd.act <- factor (LB_400.H2C5$ppfd.act, levels = c("200", "403", "850","1700" ))

plot.1.LB_400.H2C5 <- ggscatter (LB_400.H2C5 , x = "sp.fo" , y = "phiPS2", 
                      title = "LB_400.H2C5",
                      ylim=c(0,1),point.size=2,
                      cor.coef =TRUE,
                      cor.coeff.args = list(method = "pearson"),
                      facet.by = "ppfd.act",
                      fill = "ppfd.act", color = "ppfd.act",
                      add = "reg.line",
                      palette = c("orange", "red", "blue", "navyblue"))

ggexport (plot.1.LB_400.H2C5, filename = "./Figures/E1.plot.1.LB_400.H2C5.tiff",
          width = 700, height = 500)




plot.2.LB_400.H2C5 <- ggscatter (LB_400.H2C5 , x = "qp.fo" , y = "phiPS2", 
                      title = "LB_400.H2C5",
                      ylim=c(0,1),point.size=2,
                      cor.coef =TRUE,
                      cor.coeff.args = list(method = "pearson"),
                      facet.by = "ppfd.act",
                      fill = "ppfd.act", color = "ppfd.act",
                      add = "reg.line",
                       palette = c("orange", "red", "blue", "navyblue"))


ggexport (plot.2.LB_400.H2C5, filename = "./Figures/E1.plot.2.LB_400.H2C5.tiff",
          width = 700, height = 500)

### 
LB_800 <- E.1.1 %>%
  dplyr::filter (sampling == "LB_800") 

LB_800.H1C3 <- LB_800  %>%
  dplyr::filter ( genotype == "HHK1C3")

unique (LB_800.H1C3$ppfd.act )

LB_800.H1C3$ppfd.act <- factor (LB_800.H1C3$ppfd.act, levels = c("200", "403", "850","1700" ))

plot.1.LB_800.H1C3 <- ggscatter (LB_800.H1C3 , x = "sp.fo" , y = "phiPS2", 
                                 title = "LB_800.H1C3",
                                 ylim=c(0,1),point.size=2,
                                 cor.coef =TRUE,
                                 cor.coeff.args = list(method = "pearson"),
                                 facet.by = "ppfd.act",
                                 fill = "ppfd.act", color = "ppfd.act",
                                 add = "reg.line",
                                 palette = c("orange", "red", "blue", "navyblue"))


ggexport (plot.1.LB_800.H1C3, filename = "./Figures/E1.plot.1.LB_800.H1C3.tiff",
          width = 700, height = 500)



plot.2.LB_800.H1C3 <- ggscatter (LB_800.H1C3 , x = "qp.fo" , y = "phiPS2", 
                                 title = "LB_800.H1C3",
                                 ylim=c(0,1),
                                 point.size=2,
                                 cor.coef =TRUE,
                                 cor.coeff.args = list(method = "pearson"),
                                 facet.by = "ppfd.act",
                                 fill = "ppfd.act", color = "ppfd.act",
                                 add = "reg.line",
                                 palette = c("orange", "red", "blue", "navyblue"))


ggexport (plot.2.LB_800.H1C3, filename = "./Figures/E1.plot.2.LB_800.H1C3.tiff",
          width = 700, height = 500)


LB_800.H2C5 <- LB_800  %>%
  dplyr::filter ( genotype == "HHK2C5")

unique (LB_800.H2C5$ppfd.act )

LB_800.H2C5$ppfd.act <- factor (LB_800.H2C5$ppfd.act, levels = c("200", "403", "850","1700" ))


plot.1.LB_800.H2C5 <- ggscatter (LB_800.H2C5 , x = "sp.fo" , y = "phiPS2", 
                                 title = "LB_800.H2C5",
                                 ylim=c(0,1),point.size=2,
                                 cor.coef =TRUE,
                                 cor.coeff.args = list(method = "pearson"),
                                 facet.by = "ppfd.act",
                                 fill = "ppfd.act", color = "ppfd.act",
                                 add = "reg.line",
                                 palette = c("orange", "red", "blue", "navyblue"))

ggexport (plot.1.LB_800.H2C5, filename = "./Figures/E1.plot.1.LB_800.H2C5.tiff",
          width = 700, height = 500)




plot.2.LB_800.H2C5 <- ggscatter (LB_800.H2C5 , x = "qp.fo" , y = "phiPS2", 
                                 title = "LB_800.H2C5",
                                 ylim=c(0,1),point.size=2,
                                 cor.coef =TRUE,
                                 cor.coeff.args = list(method = "pearson"),
                                 facet.by = "ppfd.act",
                                 fill = "ppfd.act", color = "ppfd.act",
                                 add = "reg.line",
                                 palette = c("orange", "red", "blue", "navyblue"))


ggexport (plot.2.LB_800.H2C5, filename = "./Figures/E1.plot.2.LB_800.H2C5.tiff",
          width = 700, height = 500)

































ggline (Inv, x = "time.min", y = "phiPS2", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))

ggline (Inv, x = "time.min", y = "phiNPQ", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))


ggline (Inv, x = "time.min", y = "phiNO", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))



LB_400 <- E.1.1 %>%
          dplyr::filter (sampling == "LB_400") 

ggline (LB_400, x = "time.min", y = "phiPS2", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))


ggline (LB_400, x = "time.min", y = "phiNPQ", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))


ggline (LB_400, x = "time.min", y = "phiNO", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))

LB_800 <- E.1.1 %>%
          dplyr::filter (sampling == "LB_800") 

ggline (LB_800, x = "time.min", y = "phiPS2", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))

ggline (LB_800, x = "time.min", y = "phiNPQ", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))


ggline (LB_800, x = "time.min", y = "phiNO", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "ppfd.act",
        desc_stat = "mean_sd",
        fill = "genotype", color = "genotype",
        add = c("mean_sd"), palette = c("navyblue", "darkorange"))

#### 
unique (E.1.1$genotype)

HHK2C5  <- E.1.1 %>%
           dplyr::filter (genotype == "HHK2C5") 

HHK1C3 <- E.1.1 %>%
  dplyr::filter (genotype == "HHK1C3") 

p1 <- ggline (HHK2C5 , x = "time.min", y = "phiPS2", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "sampling",
        desc_stat = "mean_sd",
        fill = "ppfd.act", color ="ppfd.act",
        add = c("mean_sd"), palette = c("green" , " gray48", "orange", "darkred"))


p2 <- ggline (HHK1C3 , x = "time.min", y = "phiPS2", 
        ylim=c(0,1),
        point.size=2,
        facet.by = "sampling",
        desc_stat = "mean_sd",
        fill = "ppfd.act", color ="ppfd.act",
        add = c("mean_sd"), palette = c("green" , " gray48", "orange", "darkred"))

ggarrange(p1, p2, ncol = 1,labels="auto",
          nrow = 2, common.legend = TRUE)



p3 <- ggline (HHK2C5 , x = "time.min", y = "phiNPQ", 
              ylim=c(0,1),
              point.size=2,
              facet.by = "sampling",
              desc_stat = "mean_sd",
              fill = "ppfd.act", color ="ppfd.act",
              add = c("mean_sd"), palette = c("green" , " gray48", "orange", "darkred"))


p4 <- ggline (HHK1C3 , x = "time.min", y = "phiNPQ", 
              ylim=c(0,1),
              point.size=2,
              facet.by = "sampling",
              desc_stat = "mean_sd",
              fill = "ppfd.act", color ="ppfd.act",
              add = c("mean_sd"), palette = c("green" , " gray48", "orange", "darkred"))



ggarrange(p3, p4, ncol = 1,labels="auto" ,
          nrow = 2, common.legend = TRUE)


p5 <- ggline (HHK2C5 , x = "time.min", y = "phiNO", 
              ylim=c(0,1),
              point.size=2,
              facet.by = "sampling",
              desc_stat = "mean_sd",
              fill = "ppfd.act", color ="ppfd.act",
              add = c("mean_sd"), palette = c("green" , " gray48", "orange", "darkred"))


p6 <- ggline (HHK1C3 , x = "time.min", y = "phiNO", 
              ylim=c(0,1),
              point.size=2,
              facet.by = "sampling",
              desc_stat = "mean_sd",
              fill = "ppfd.act", color ="ppfd.act",
              add = c("mean_sd"), palette = c("green" , " gray48", "orange", "darkred"))



ggarrange(p5, p6, ncol = 1,labels="auto" ,
          nrow = 2, common.legend = TRUE)











ggerrorplot (Inv , x = "time.min", y = "phiPS2", 
             facet.by = "genotype" ,
             ylim=c(0,1),
             fill = "ppfd.act", color = "ppfd.act",
             desc_stat = "mean_sd",
             error.plot = "errorbar",            # Change error plot type
             add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red")) 



HHK2C5.5.5 <- E.1.1 %>%
          dplyr::filter (genotype == "HHK2C5") %>%
          dplyr::filter (time.min == 5.5)

##

ggerrorplot (HHK2C5.5.5 , x = "ppfd.act", y = "phiPS2", 
             #facet.by = "ppfd.act" ,
             ylim=c(0,1),
             fill = "sampling", color = "sampling",
             desc_stat = "mean_sd",
             error.plot = "errorbar",            # Change error plot type
             add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "red"))


ggline (HHK2C5.5.5, x = "ppfd.act", y = "phiPS2", 
       ylim=c(0,1),point.size=2,
       desc_stat = "mean_sd",
       fill = "sampling", color = "sampling",
       add = c("mean_sd"), palette = c("navyblue", "gray48", "red"))



ggline (HHK2C5.5.5, x = "ppfd.act", y = "phiPS2", 
        ylim=c(0,1),point.size=2,
        desc_stat = "mean_sd",
        fill = "sampling", color = "sampling",
        add = c("mean_sd"), palette = c("navyblue", "gray48", "red"))











ggline (HHK2C5.5.5, x = "dose", y = "len", color = "supp",
        



ggerrorplot (GT_arcadida , x = "time.min", y = "qp.fo", 
             facet.by = "actinica",
             ylim=c(0,1),
             fill = "QIc", color = "QIc",
             desc_stat = "mean_sd",
             #error.plot = "errorbar",            # Change error plot type
             #add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red")) 


ggerrorplot (GT_arcadida , x = "time.min", y = "sp.fo", 
             facet.by = "actinica",
             ylim=c(0,1),
             fill = "QIc", color = "QIc",
             desc_stat = "mean_sd",
             #error.plot = "errorbar",            # Change error plot type
             #add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red")) 


ggerrorplot (GT_arcadida , x = "time.min", y = "phiNPQ", 
             facet.by = "actinica",
             ylim=c(0,1),
             fill = "QIc", color = "QIc",
             desc_stat = "mean_sd",
             #error.plot = "errorbar",            # Change error plot type
             #add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red")) 


ggerrorplot (GT_arcadida , x = "time.min", y = "npq1.fo", 
             facet.by = "actinica",
             ylim=c(0,1),
             fill = "QIc", color = "QIc",
             desc_stat = "mean_sd",
             #error.plot = "errorbar",            # Change error plot type
             #add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red")) 


ggerrorplot (GT_arcadida , x = "time.min", y = "npq2.fo", 
             facet.by = "actinica",
             ylim=c(0,1),
             fill = "QIc", color = "QIc",
             desc_stat = "mean_sd",
             #error.plot = "errorbar",            # Change error plot type
             #add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red")) 



ggerrorplot (GT_arcadida , x = "time.min", y = "phiNO", 
             facet.by = "actinica",
             ylim=c(0,1),
             fill = "QIc", color = "QIc",
             desc_stat = "mean_sd",
             #error.plot = "errorbar",            # Change error plot type
             #add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red"))


ggerrorplot (GT_arcadida , x = "time.min", y = "no1.fo", 
             facet.by = "actinica",
             ylim=c(0,1),
             fill = "QIc", color = "QIc",
             desc_stat = "mean_sd",
             #error.plot = "errorbar",            # Change error plot type
             #add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red"))



ggerrorplot (GT_arcadida , x = "time.min", y = "no2.fo", 
             facet.by = "actinica",
             ylim=c(0,1),
             fill = "QIc", color = "QIc",
             desc_stat = "mean_sd",
             #error.plot = "errorbar",            # Change error plot type
             #add = "mean",
             #position = position_dodge(1) ,
             palette = c("navyblue", "gray48", "black", "red"))



       add = c("mean_se", "jitter"), palette = c("#00AFBB", "#E7B800"))


GT.1_relax.S3 <- run_relax_analisis ( dt = PAM.param.s3.a, sampling ="GT.1",protocol = "s3", rep ="si")


GT.1_S3.1 <- run_relax_analisis ( dt = PAM.param.s3.1.a, sampling ="GT.1",protocol = "s3.1", rep ="si")


GT.1_S3.2 <- run_relax_analisis ( dt = PAM.param.s3.2.a, sampling ="GT.1",protocol = "s3.2", rep ="si")




preliminar_S3.1 <- run_quenching_analisis ( dt = PAM.param.1, sampling ="preliminar",protocol = "Simondi3.1")

preliminar_S3.2 <- run_quenching_analisis ( dt = PAM.param.1, sampling ="preliminar",protocol = "Simondi3.2")


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
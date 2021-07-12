######################################################################################
# Codigo para el analisis de la base de datos de Fluorescencia en Cebada           #
#                                                                                    #
# Gaston Quero - Sebastian Simondi                                                   #          
# 18/3/2021                                                                          #
######################################################################################

#####
# los datos correspoden a .......

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
library("Hmisc")
library("PerformanceAnalytics")
library (GGally)
library(ggcorrplot)
library (readr)
library (igraph)
library (ggraph)

#install.packages("PerformanceAnalytics")
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


# reingreso los datos y saco el ultimo ciclo #### 


Phi.PSII.SP.1b <- read_delim (file ="./Data/procdata/Phi.PSII.SP.1b.txt", 
                  delim =";", quote = "\"", escape_backslash = FALSE,
                  escape_double = TRUE, col_names = TRUE, col_types = NULL,
                  na = "NA")

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

Phi.PSII.SP.1.ciclo.q <- read_delim (file ="./Data/procdata/Phi.PSII.SP.1.ciclo.q.txt", 
                              delim =";", quote = "\"", escape_backslash = FALSE,
                              escape_double = TRUE, col_names = TRUE, col_types = NULL,
                              na = "NA")



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
  
unique ( Phi.PSII.SP.1.ciclo.q.1$amb.ppfd  )

LB800.1700 <- Phi.PSII.SP.1.ciclo.q.1 %>%
  dplyr::filter (amb.ppfd ==  "LB800") %>%
  dplyr::filter ( ppfd.act  == 1700)

potX <- unique (LB800.1700$filt.pot)

Phi.PSII.SP.1.ciclo.q.1 <-  Phi.PSII.SP.1.ciclo.q.1 %>%
  dplyr::filter (!filt.pot %in% potX)


# primer grafo con la particion principal de energia 

phis.1 <- c ("phiPS2", "phiNPQ", "phiNO" ) 

# selecciono las variables a correlacionar
phis.cor <- Phi.PSII.SP.1.ciclo.q.1 %>%
             dplyr::select (all_of(phis.1))

# hago las correlaciones 
res.phi.2 <- rcorr (as.matrix (phis.cor))
res.phi.2 

# esta es la funcion que transforma la matriz de correlaciones en un data.frame
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# genero el data frame con las correlaciones
res.phi.3 <- flattenCorrMatrix (res.phi.2$r, res.phi.2$P) ### esto lo puedo usar para el grafo


# To visualize the network graph, we need to create two data frames from the data sets:
  
# nodes list: containing nodes labels and other nodes attributes
# edges list: containing the relationship between the nodes. It consists of the edge list and any additional edge attributes.


#Create nodes list

#Create nodes list

res.phi.3 <- as_tibble (res.phi.3)
#  Get distinct source names
origen <- res.phi.3 %>%
          dplyr::distinct (row) %>%
          dplyr::rename(label = row)

# Get distinct destination names
destino <- res.phi.3 %>%
           dplyr::distinct (column) %>%
           dplyr::rename (label = column)


# Join the two data to create node
# Add unique ID for each country
vertices <- full_join (origen, destino, by = "label") 

vertices <- vertices %>%
            dplyr::mutate(id = 1:nrow(vertices)) %>%
            select(id, everything())

head (vertices , 3)

#Take the phone.call data, which are already in edges list format,
# showing the connection between nodes. 
# Rename the column “n.call” to “weight”.
# Join the node IDs to the edges list data
# Do this for the “source” column and rename the id column that are brought over from nodes. New name: “from”.
# Do this for the “destination” column and rename the id column. New name: “to”
# Select only the columns “from” and “to” in the edge data. 
# We don’t need to keep the column “source” and “destination” containing the names of countries. 
# These information are already present in the node data.

# Rename the n.call column to weight
res.phi.3 <- res.phi.3 %>%
              rename(weight = cor)

# (a) Join nodes id for source column
aristas <- res.phi.3 %>% 
           left_join (vertices, by = c("row" = "label")) %>% 
           dplyr::rename (from = id) 

# (b) Join nodes id for destination column
#fujoE <- tibble (energia = c("F.C", "F.Cb", "C.Cb" ) )

aristas <- aristas %>% 
         left_join (vertices, by = c("column" = "label")) %>% 
         dplyr::rename(to = id) %>%
         dplyr::mutate (weight = round (weight, 2))%>%
         dplyr::select (from, to, weight)




head(aristas, 3)

set.seed(123)

#### igraph #######
## Create an igraph network object:
# Key R function: graph_from_data_frame().

# Key arguments:
# d: edge list
# vertices: node list
# directed: can be either TRUE or FALSE depending on whether the data is directed or undirected.

net.principal.phi <- graph_from_data_frame(
  d = aristas, vertices = vertices, 
  directed = FALSE
)

# Create a network graph with igraph

set.seed(321)
plot (net.principal.phi, 
      edge.arrow.size = 0.2,
     layout = layout_with_graphopt)



# tidygraph and ggraph

#Create a network object using tidygraph:
#  Key function: tbl_graph().
# key arguments: nodes, edges and directed.
library(tidygraph)

net.tidy.principal.phi <- tbl_graph (
  nodes = vertices, edges = aristas, directed = FALSE
)

# Not specifying the layout - defaults to "auto"


layout.x <- create_layout ( net.tidy.principal.phi, layout = 'auto')

layout.x$x[1] <- 0
layout.x$x[2] <- -0.86
layout.x$x[3] <- 0.86

layout.x$y[1] <- 0.5
layout.x$y[2] <- 0
layout.x$y[3] <- 0


ggraph (layout.x) + 
     geom_edge_link (alpha = 0.25,
                       aes(width = abs( weight), label = weight, repel = TRUE)) +
  geom_node_point(aes ( size =3, alpha =0.5)) +
  geom_node_text(aes(label = label), repel = TRUE) +
 theme_graph() 

attributes(layout.x)



# primer grafo
#phis.1 <- c ("phiPS2", "phiNPQ", "phiNO" ) 

phis.2 <- c ( "qp.fo" , "qs.fo" , "phiNO.psII" ,
              "phiNO.basal" , "phiNPQfast", "phiNPQslow" ) 

# selecciono las variables a correlacionar
phis.cor <- Phi.PSII.SP.1.ciclo.q.1 %>%
  dplyr::select (all_of(phis.2))

# hago las correlaciones 
res.phi.2 <- rcorr (as.matrix (phis.cor))
res.phi.2 

# esta es la funcion que transforma la matriz de correlaciones en un data.frame
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# genero el data frame con las correlaciones
res.phi.3 <- flattenCorrMatrix (res.phi.2$r, res.phi.2$P) ### esto lo puedo usar para el grafo


# To visualize the network graph, we need to create two data frames from the data sets:

# nodes list: containing nodes labels and other nodes attributes
# edges list: containing the relationship between the nodes. It consists of the edge list and any additional edge attributes.


#Create nodes list

#Create nodes list

res.phi.3 <- as_tibble (res.phi.3)
#  Get distinct source names
origen <- res.phi.3 %>%
  dplyr::distinct (row) %>%
  dplyr::rename(label = row)

# Get distinct destination names
destino <- res.phi.3 %>%
  dplyr::distinct (column) %>%
  dplyr::rename (label = column)


# Join the two data to create node
# Add unique ID for each country
vertices <- full_join (origen, destino, by = "label") 

vertices <- vertices %>%
  dplyr::mutate(id = 1:nrow(vertices)) %>%
  select(id, everything())

head (vertices , 3)

#Take the phone.call data, which are already in edges list format,
# showing the connection between nodes. 
# Rename the column “n.call” to “weight”.
# Join the node IDs to the edges list data
# Do this for the “source” column and rename the id column that are brought over from nodes. New name: “from”.
# Do this for the “destination” column and rename the id column. New name: “to”
# Select only the columns “from” and “to” in the edge data. 
# We don’t need to keep the column “source” and “destination” containing the names of countries. 
# These information are already present in the node data.

# Rename the n.call column to weight
res.phi.3 <- res.phi.3 %>%
  rename(weight = cor)

# (a) Join nodes id for source column
aristas <- res.phi.3 %>% 
  left_join (vertices, by = c("row" = "label")) %>% 
  dplyr::rename (from = id) 

# (b) Join nodes id for destination column
#fujoE <- tibble (energia = c("F.C", "F.Cb", "C.Cb" ) )

aristas <- aristas %>% 
  left_join (vertices, by = c("column" = "label")) %>% 
  dplyr::rename(to = id) %>%
  dplyr::mutate (weight = round (weight, 2))%>%
  dplyr::select (from, to, weight)




head(aristas, 3)

set.seed(123)

#### igraph #######
## Create an igraph network object:
# Key R function: graph_from_data_frame().

# Key arguments:
# d: edge list
# vertices: node list
# directed: can be either TRUE or FALSE depending on whether the data is directed or undirected.

net.principal.phi <- graph_from_data_frame(
  d = aristas, vertices = vertices, 
  directed = FALSE
)

# Create a network graph with igraph

set.seed(321)
plot (net.principal.phi, 
      edge.arrow.size = 0.2,
      layout = layout_with_graphopt)



# tidygraph and ggraph

#Create a network object using tidygraph:
#  Key function: tbl_graph().
# key arguments: nodes, edges and directed.
#library(tidygraph)

net.tidy.principal.phi <- tbl_graph (
  nodes = vertices, edges = aristas, directed = FALSE
)

# Not specifying the layout - defaults to "auto"



layout.x <- create_layout (net.tidy.principal.phi, layout = 'auto')

layout.x$x[1] <- 2
layout.x$x[2] <- 4
layout.x$x[3] <- 1
layout.x$x[4] <- 5
layout.x$x[5] <- 2
layout.x$x[6] <- 4




layout.x$y[1] <- 3
layout.x$y[2] <- 3
layout.x$y[3] <- 2
layout.x$y[4] <- 2
layout.x$y[5] <- 1
layout.x$y[6] <- 1

ggraph (layout.x) + 
  geom_edge_link ( alpha = 0.25, 
                   aes(width = abs( weight), label = weight, repel = TRUE)) +
  geom_node_point(aes ( size =3, alpha =0.5)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  theme_graph() 

attributes(layout.x)







create_layout.myclass()







set.seed(123)
ggraph(net.tidy.principal.phi , layout = "lgl") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.5) + 
  #scale_edge_width (range = c(0, 1)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "correlation") +
  theme_graph()

set.seed(123)
ggraph(net.tidy.principal.phi, layout = "linear") + 
  geom_edge_arc(aes(width = weight), alpha = 0.8) + 
  #scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "corr") +
  theme_graph()+
  theme(legend.position = "top")




# Coord diagram, circular
ggraph(net.tidy.principal.phi, layout = "linear", circular = TRUE) + 
  geom_edge_arc(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "corr") +
  theme_graph()+
  theme(legend.position = "top")







#star', 
# 'circle', 
#'gem', 'dh', 
#''graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl'.


nodes <- nam

#%>%
 #        dplyr::mutate ( id = row.names(nam)) %>%
  #       dplyr::select (id , name)


ties <- res.phi.3 %>%
        dplyr::mutate (from = row) %>%
        dplyr::mutate (to= column) %>%
        dplyr::mutate (weight = cor) %>%
        dplyr::select (from, to, weight) %>%
        dplyr::mutate (weight = round (weight, 3))


ties [ties == "phiPS2"] <- 1
#ties [ties == "qp.fo"]  <- 2
#ties [ties == "qs.fo"]  <- 3
ties [ties == "phiNPQ"]  <- 2
ties [ties == "phiNO"]  <- 3
#ties [ties == "phiNO.psII"]  <- 6
#ties [ties == "phiNO.basal"]  <- 7
#ties [ties == "phiNPQfast"]  <- 8
#ties [ties == "phiNPQslow"]  <- 9
# Print nodes
print (nodes)

# Print ties
print (ties)


 # Make the network from the data frame ties and print it
g <- graph_from_data_frame (ties, directed = TRUE, vertices = nodes)
g

# Give the name "Madrid network" to the network 
g$name <- "SP.1"
g$name


# Add node attribute id and print the node `id` attribute
V(g)$id <- c ("phiPS2", "phiNPQ", "phiNO" )
V(g)$id

# Print the tie `weight` attribute
E(g)$weight

plot (g)

ggraph (g, layout = 'fr') +
geom_node_point(color = "navyblue", size = 5) +
    geom_edge_link ( alpha = 0.25, 
                     aes(width = abs( weight), label = weight,
                                       start_cap = label_rect(node1.name),
                                       end_cap = label_rect(node2.name)))+
  geom_node_point(color = "navyblue", size = 5) + 
  geom_node_text(aes(label = name),  repel = TRUE)+
    coord_fixed() +
  theme_graph()













ggraph(g, layout = 'graphopt') + 
  geom_edge_link(aes(start_cap = label_rect(node1.name),
                     end_cap = label_rect(node2.name)), 
                 arrow = arrow(length = unit(4, 'mm'))) + 
  geom_node_text(aes(label = name))

ggraph(g, layout = 'graphopt') + 
  geom_edge_link(aes(label = weight), 
                 label_dodge = unit(2.5, 'mm'),
                 arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm')) + 
  geom_node_point(size = 5) +
  geom_node_text (aes(label = id), repel = TRUE)+
  coord_fixed()
geom_edge_link(aes(label =  weight), 
       
ggraph(g, layout = 'graphopt') + 
            angle_calc = 'along',
                 label_dodge = unit(2.5, 'mm'),
                 arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm')) + 
  geom_node_point(size = 5) +
geom_node_text (aes(label = id), repel = TRUE)+
  coord_fixed()


cor.pos <- degree(delete_edges(g, which (E(g)$weight > 0 )), 
                  mode = 'in')

cor.neg <- degree(delete_edges(g, which (E(g)$weight < 0 )), 
                  mode = 'in')


ggraph(g, layout = 'hive') + 
  geom_edge_hive(aes(colour = weight)) + 
  coord_fixed()


# Explore the set of nodes
V(g)

# Print the number of nodes
vcount (g)

# Explore the set of ties
E(g)

# Print the number of ties
ecount (g)

# Give the name "Madrid network" to the network 
g$name <- "SP.1"
g$name


# Add node attribute id and print the node `id` attribute
V(g)$id <- c ("phiPS2", "phiNPQ", "phiNO" )
V(g)$id

# Print the tie `weight` attribute
E(g)$weight

# Print the network and spot the attributes
g


plot (g)
# Visualize the network with the Kamada-Kawai layout 
ggraph(g, layout = "auto") + 
  # Add an edge link geometry mapping transparency to weight 
  geom_edge_link (aes(colour = weight)) + 
  # Add a node point geometry
  geom_node_point(aes (size = 1) )+
  geom_node_text (aes(label = id), repel = TRUE)


# Visualize the network in a circular layout
ggraph(g, layout = "in_circle", circular = FALSE ) + 
  # Map tie transparency to its weight
  geom_edge_link(aes(alpha = weight)) + 
  geom_edge_link (aes(size = weight)) +
  geom_node_point()+
  geom_node_text (aes(label = id), repel = TRUE)


# Change the layout so points are on a grid
ggraph(g, layout = "on_grid", circular = FALSE ) + 
  # Map tie transparency to its weight
  geom_edge_link(aes(alpha = weight)) + 
  geom_node_point()+
  geom_node_text (aes(label = id), repel = TRUE)


ggraph(g) + 
  geom_edge_link(aes(colour = weight, alpha = weight))
  
  geom_edge_link(aes(colour = factor(weight))) + 
  geom_node_point()


ggraph(g, layout = 'kk') + 
  geom_edge_link(aes(colour = factor(weight))) + 
  geom_node_point()


ggraph (hairball, layout = 'hive', axis = 'pop_devel', sort.by = 'popularity') + 
  geom_edge_hive(aes(colour = year)) + 
  geom_axis_hive() + 
  coord_fixed()

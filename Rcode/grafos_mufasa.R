if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/navdata")

install.packages(
  c("tidyverse", "igraph", "tidygraph", "ggraph")
)
library("navdata")
data("phone.call")
head(phone.call, 3)

#Create nodes list
#  Get distinct source names
sources <- phone.call %>%
           distinct(source) %>%
           rename(label = source)

# Get distinct destination names
destinations <- phone.call %>%
                distinct(destination) %>%
                rename(label = destination)


# Join the two data to create node
# Add unique ID for each country
nodes <- full_join (sources, destinations, by = "label") 
nodes <- nodes %>%
         mutate(id = 1:nrow(nodes)) %>%
         select(id, everything())


head(nodes, 3)

# Create edges list
# Key steps:
  
# Take the phone.call data, which are already in edges list format,
# showing the connection between nodes. 
# Rename the column “n.call” to “weight”.
# Join the node IDs to the edges list data
# Do this for the “source” column and rename the id column that are brought over from nodes. New name: “from”.
# Do this for the “destination” column and rename the id column. New name: “to”
# Select only the columns “from” and “to” in the edge data. 
# We don’t need to keep the column “source” and “destination” containing the names of countries. 
# These information are already present in the node data.

# Rename the n.call column to weight
phone.call <- phone.call %>%
              rename(weight = n.call)

# (a) Join nodes id for source column
edges <- phone.call %>% 
         left_join(nodes, by = c("source" = "label")) %>% 
         rename(from = id)

# (b) Join nodes id for destination column
edges <- edges %>% 
         left_join (nodes, by = c("destination" = "label")) %>% 
         rename(to = id)

# (c) Select/keep only the columns from and to
edges <- select (edges, from, to, weight)
head(edges, 3)


#### igraph #######
## Create an igraph network object:
# Key R function: graph_from_data_frame().

# Key arguments:
# d: edge list
# vertices: node list
# directed: can be either TRUE or FALSE depending on whether the data is directed or undirected.


net.igraph <- graph_from_data_frame(
  d = edges, vertices = nodes, 
  directed = TRUE
)

# Create a network graph with igraph
set.seed(123)
plot (net.igraph, edge.arrow.size = 0.2,
      layout = layout_with_graphopt)

?plot.igraph

# tidygraph and ggraph

#Create a network object using tidygraph:
#  Key function: tbl_graph().
# key arguments: nodes, edges and directed.
library(tidygraph)
net.tidy <- tbl_graph(
  nodes = nodes, edges = edges, directed = TRUE
)

# Visualize network using ggraph
# Key functions:
  
#geom_node_point(): Draws node points.

#geom_edge_link(): Draws edge links. 

########################################################################
#To control the width of edge line according to the weight variable,
# specify the option aes(width = weight), where, the weight specify the number of phone.call sent along each route. 
# In this case, you can control the maximum and minimum width of the edges, by using the function scale_edge_width() to set the range (minimum and maximum width value).
# For example: scale_edge_width(range = c(0.2, 2)).
#####################################################################

# geom_node_text(): Adds text labels for nodes, by specifying the argument aes(label = label). 
# To avoid text overlapping, indicate the option repel = TRUE.

# labs(): Change main titles, axis labels and legend titles.

# Create a classic node-edge diagrams. 
# Possible values for the argument layout include:
 # 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl'.

set.seed(123)
ggraph(net.tidy, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "phone.call") +
  theme_graph()


# Arc diagram
ggraph(net.tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Number of calls") +
  theme_graph()+
  theme(legend.position = "top") 


##### Thomas Lin Pedersen ###3

library(ggraph)
library(igraph)
graph <- graph_from_data_frame(highschool)

# Not specifying the layout - defaults to "auto"
ggraph(graph) + 
  geom_edge_link(aes(colour = factor(year))) + 
  geom_node_point()

layout <- create_layout(graph, layout = 'drl')
ggraph(layout) + 
  geom_edge_link(aes(colour = factor(year))) + 
  geom_node_point()
attributes(layout)

gr <- as_tbl_graph(highschool)

ggraph(gr, layout = 'fabric') +
  geom_node_range()


source("global.R")
source("functions.R")
#source("prepare.R")

#### 0. SUBSET PATENTS TO GIVEN USPC CLASS AND APPLICATION YEAR (OPTIONAL) ####

load("data/invpat.RData")

## keep patents with the USPC class 438: Semiconductor device manufacturing: process 
## and application year == 2000
invpat <- invpat_raw[!is.na(invpat_raw$Class) & invpat_raw$Class=="438" &
                       !is.na(invpat_raw$AppYearStr) & invpat_raw$AppYearStr==2000,]
rm(invpat_raw)
gc()

save(invpat,file="invpat_sub.RData")

#### 1. COLLABORATION NETWORK (LAYER 1) CONSTRUCTION ####

load("data/invpat_sub.RData")
invpat <- invpat[,c("Patent","Invcode","LocCode","Lat","Lon")]

## NB: the input data containt the following columns:
##     1) unique identifier of a node (inventor);
##     2) unique identifier of a collaboration team (patent);
##     3) unique identifier of a location (pair of longitude and latitude)
##     4) longitude and latitude of an inventor

## for each team (patent) create a list of all possible unordered pairs of its members (coauthors)
pat_unique <- unique(invpat$Patent)
edgelist <- vector()

p <- 1
for(i in 1:length(pat_unique)){
  subset <- invpat$Invcode[invpat$Patent == pat_unique[i]]
  if(length(subset)>1){edgelist <- append(edgelist,c(combn(subset,2)))}
  
  progress(p,length(c(1:length(pat_unique))))
  Sys.sleep(0.01)
  if (p == length(c(1:length(pat_unique)))) cat("\nList of collaboration pairs is created!")
  p <- p+1
}

## create a network of inventors (igraph object) from a list of collaboration pairs
network <- igraph::graph(edgelist,directed=FALSE)
network <- igraph::simplify(network)
rm(pat_unique,edgelist,subset,i)

## add location coordinates as attributes of nodes
V(network)$Lat <- invpat$Lat[match(V(network)$name,invpat$Invcode)]
V(network)$Lon <- invpat$Lon[match(V(network)$name,invpat$Invcode)]

dir.create("./aux")
save(network,file="aux/collaboration_network.RData")
  
#### 2. GEOGRAPHICAL NETWORK (LAYER 2) CONSTRUCTION ####

load("data/invpat_sub.RData")
load("aux/collaboration_network.RData")
invpat <- invpat[,c("Patent","Invcode","LocCode","Lat","Lon")]

## NB: as opposed to the collaboration network, where each node is an inventor, the geographical
##     network is constructed among unique geographical locations (for simplification)

## NB: each inventor (node) is supposed to have a unique location in the dataset

## take a set of unique locations
loc <- unique(invpat[,c("LocCode","Lon","Lat")])
loc <- loc[order(loc$LocCode),]

## construct a matrix of distances (calculated using the Haversive formula, in km) between any two unique locations
dist_matrix <- unname(apply(loc[,c("Lon","Lat")], 1, FUN=function(X) geosphere::distHaversine(X, loc[,c("Lon","Lat")])/1000))

## loop over a range of distance thresholds (maximum distance separating the node from its closest geographical neighbor)
dist_threshold <- c(seq(1,50,1),seq(55,100,5),seq(110,200,10))

dir.create("./geocomponents")
p <- 1
for(d in dist_threshold){
  
  ## convert distance matrix into adjacency matrix (based on the distance threshold)
  ## convert adjacency matrix into the network (igraph object)
  geo_network <- igraph::graph_from_adjacency_matrix((dist_matrix<=d & dist_matrix>0)*1, mode="upper")
  
  ## add location identifier as an attribute of nodes
  V(geo_network)$name <- loc$LocCode
  V(geo_network)$LocCode <- loc$LocCode
  
  ## use the igraph function "clusters" to split the geographical network into geo-components (group the locations)
  geo_components <- igraph::clusters(geo_network)
  geo_components <- data.frame(LocCode = V(geo_network)$LocCode,
                               GeoComponentCode = geo_components$membership)
  
  ## match locations in the geo-components with the inventors
  inv_geo_component <- merge(unique(invpat[,c("Invcode","LocCode")]),geo_components)
  
  ## keep only those inventors who appear on a collaboration network (i.e. have at least one collaboration with others)
  inv_geo_component <- inv_geo_component[inv_geo_component$Invcode %in% V(network)$name,]
  
  ## create a list (over geo-components) of vectors with inventor codes belonging to separate geo-components
  inv_geo_component <- inv_geo_component[order(inv_geo_component$GeoComponentCode),]
  inv_geo_component <- tapply(inv_geo_component$Invcode, inv_geo_component$GeoComponentCode, function(x) x )
  
  ## calculate the size (number of members) of each geo-component
  size_geo_component <- unlist(lapply(inv_geo_component, function(x) length(x) ))
  
  save(list=c("inv_geo_component","size_geo_component","d"),file=paste0("geocomponents/geo_components_",d,".RData"))
  
  progress(p,length(dist_threshold))
  Sys.sleep(0.01)
  if (p == length(dist_threshold)) cat("\nSet of geographical networks (and their components) is created!")
  p <- p+1
}

#### 3. IDENTIFY CLUSTERS ####
rm(list=ls())
source("functions.R")

load("aux/collaboration_network.RData")

## loop over a range of distance thresholds (maximum distance separating the node from its closest geographical neighbor)
dist_threshold <- c(seq(1,50,1),seq(55,100,5),seq(110,200,10)) # threshold km for identifying geographic neighborhood

## loop over a range of size thresholds (minimum number of cluster members)
size_threshold <- c(seq(5,50,5),seq(60,200,10))

dir.create("./clusters")
for(d in dist_threshold){
  
  load(paste0("geocomponents/geo_components_",d,".RData"))
  
  for(s in size_threshold){
    inv_clusters <- list()
    inx <- vector()

    ## consider geo-components of the geographical network that are larger than threshold
    inv_geo_component_large <- inv_geo_component[size_geo_component>s]
    
    if(length(inv_geo_component_large)>0){
      
      ## loop over the considered geo-components
      p <- 1
      for(i in 1:length(inv_geo_component_large)){
        
        ## code 0 (iteration 0)
        sub_network_0 <- igraph::induced_subgraph(network,V(network)[inv_geo_component_large[[i]]])
        connected_components_0 <- igraph::clusters(sub_network_0)
        connected_components_0 <- unlist_component(connected_components_0)
        geo_components_x <- rm_null(lapply(connected_components_0, function(x) geo_components(sub_network_0,x,d) ))
        
        while(length(geo_components_x)>0){
          ## off
          if(length(inx_off(geo_components_x))>0){
            inv_clusters[(length(inv_clusters)+1):(length(inv_clusters)+length(inx_off(geo_components_x)))] <- 
              lapply(geo_components_x[inx_off(geo_components_x)],function(x)names(x$membership))
          }
          ## code 1 (on)
          geo_components_x <- unlist(lapply(geo_components_x[inx_on(geo_components_x)],function(x)unlist_component(x)),recursive=FALSE)
          ## code 2
          connected_components_x <- lapply(geo_components_x,function(x)connected_components(sub_network_0,x))
          ## off
          if(length(inx_off(connected_components_x))>0){
            inv_clusters[(length(inv_clusters)+1):(length(inv_clusters)+length(inx_off(connected_components_x)))] <- 
              lapply(connected_components_x[inx_off(connected_components_x)],function(x)names(x$membership))
          }
          ## code 3 (on)
          connected_components_x <- unlist(lapply(connected_components_x[inx_on(connected_components_x)],function(x)unlist_component(x)),recursive=FALSE)
          ## code 4
          geo_components_x <- lapply(connected_components_x,function(x)geo_components(sub_network_0,x,d))
        }
        
        progress(p,length(inv_geo_component_large))
        Sys.sleep(0.01)
        if (p == length(inv_geo_component_large)) progress_message(d,s,i)
        p <- p+1
      }
    }else warning(paste("There are no large enough geo-components to start with for distance threshold:",d,"and size threshold:",s))
    save(list=c("inv_clusters","dist_threshold","size_threshold"),file=paste0("clusters/clusters_",d,"_",s,".RData"))
  }
}
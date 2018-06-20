# FUNCTIONS FOR CLUSTER IDENTIFICATION ###############################################################

## this function converts connected components in 'list' format into geo-components in 'igraph' format
geo_components <- function(g,inv_connected_components,distance_threshold){
  
  ## extract inventor coordinates from the attributes of collaboration network
  connected_component <- induced.subgraph(g,V(g)[match(inv_connected_components,V(g)$name)])
  coord_matrix <- cbind(V(connected_component)$Lon,V(connected_component)$Lat)
  
  ## create a list of edges (pairs of inventors who are closer than d kilometers from each other)
  inv_vector <- V(connected_component)$name
  edgelist <- matrix(NA,0,2)
  for(z in 1:(nrow(coord_matrix)-1)){
    x <- distHaversine(coord_matrix[z,],coord_matrix[-c(1:z),])/1000
    y <- inv_vector[-c(1:z)]
    if(length(which(x<=distance_threshold))>0){
      edgelist <- rbind(edgelist,cbind(rep(inv_vector[z],length(which(x<=distance_threshold))),y[x<=distance_threshold]))
    }
  }
  
  ## if there is at least one such pair, construct a geographical network and split it into geo-components
  if(nrow(edgelist)>0){return(clusters(graph.data.frame(data.frame(edgelist),vertices=inv_vector)))}
}

## this function converts geo-components in 'list' format into connected components in 'igraph' format
connected_components <- function(g,inv_geo_components){
  geo_component <- induced.subgraph(g,V(g)[match(inv_geo_components,V(g)$name)])
  return(clusters(geo_component))
}

## this function converts the connected components / geo-components from 'igraph' format into the 'list' format
unlist_component <- function(components){
  return(lapply(which(components$csize>s),function(x)names(components$membership[components$membership==x])))
}

## this function extracts connected components that entirely coincide with geo-components, i.e. are identified as clusters
inx_off <- function(components){
  which(unlist(lapply(components,function(x)x$no))==1 & unlist(lapply(rm_null(components),function(x)max(x$csize)))>s)
}

## this function extracts connected components that can be splited into more than one geo-component, i.e. should proceeded to the next iteration
inx_on <- function(components){
  which(unlist(lapply(components,function(x)x$no))>1 & unlist(lapply(rm_null(components),function(x)max(x$csize)))>s)
}

## this function removes NULL elements of the list
rm_null <- function(components){
  components[!sapply(components, is.null)]
}

## this is an auxiliarty function used for progress reporting
progress_message <- function(d,s,i){
  cat(paste("\nClusters for distance threshold:",d,"and size threshold:",s,"are identified!\n",
            paste(rep(" ",15+nchar(i)*2),collapse="")))
}

# FUNCTIONS FOR CLUSTER VISUALIZATION ################################################################
make_circles <- function(centers, radius, nPoints = 100){
  ## centers: the data frame of centers with ID
  ## radius: radius measured in kilometer
  
  meanLat <- mean(centers$Lat)
  ## length per longitude changes with lattitude, so need correction
  radiusLon <- radius /111 / cos(meanLat/57.3) 
  radiusLat <- radius / 111
  circleDF <- data.frame(ID = rep(centers$LocCode, each = nPoints))
  angle <- seq(0,2*pi,length.out = nPoints)
  
  circleDF$lon <- unlist(lapply(centers$Lon, function(x) x + radiusLon * cos(angle)))
  circleDF$lat <- unlist(lapply(centers$Lat, function(x) x + radiusLat * sin(angle)))
  return(circleDF)
}
######################################################################################################
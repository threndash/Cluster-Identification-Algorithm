source("global.R")
source("functions.R")

#### 0. PREPARE DATA FOR VISUALIZATION ####

load("data/invpat_sub.RData")
load("aux/collaboration_network.RData")

## unique inventors with their coordinates
inv_loc <- unique(invpat[,c("Invcode","Lat","Lon")])

## keep locations within the main North American continent
inv_loc <- inv_loc[inv_loc$Lat<=50 &
                     inv_loc$Lat>=25 &
                     inv_loc$Lon<=-65 &
                     inv_loc$Lon>=-125,]
inv_loc$LocCode <- as.numeric(factor(paste(inv_loc$Lat,inv_loc$Lon)))

#### 1. MERGE WITH THE RESULTS FROM CLUSTER IDENTIFICATION ####

## set the size and distance thresholds here
d = 9
s = 10

load(paste0("clusters/clusters_",d,"_",s,".RData"))

## convert list of clusters into data.frame of inventor codes and cluster IDs
inv_clusters <- data.frame(Invcode = unlist(inv_clusters),
                           Cluster_ID = rep(seq_len(length(inv_clusters)),sapply(inv_clusters,length)))

## merge it with all inventors in the dataset
inv_loc_clusters <- merge(inv_loc,inv_clusters)

#### 2. VISUALIZE CLUSTERS ####
## create maps of locations and networks of collaborations for a given set of clusters

dir.create("./maps")
dir.create("./networks")

p <- 1
for(i in unique(inv_loc_clusters$Cluster_ID)){
  
  cluster_members <- inv_loc_clusters[inv_loc_clusters$Cluster_ID %in% i,]
  links <- induced.subgraph(network,cluster_members$Invcode)

  coords <- as.data.frame(cbind(V(network)$Lon[match(get.edgelist(links)[,1],V(network)$name)],
                                V(network)$Lat[match(get.edgelist(links)[,1],V(network)$name)],
                                V(network)$Lon[match(get.edgelist(links)[,2],V(network)$name)],
                                V(network)$Lat[match(get.edgelist(links)[,2],V(network)$name)]))
  names(coords) <- c("x1","y1","x2","y2")
  myCircles_clu <- make_circles(cluster_members[,c("Lat","Lon","LocCode")], d)
  
  x_limits <- c(min(myCircles_clu$lon)-0.1, max(myCircles_clu$lon)+0.1)
  y_limits <- c(min(myCircles_clu$lat)-0.1, max(myCircles_clu$lat)+0.1)
  
  while(!exists("background_layer")){
    background_layer <- qmap(location = c(lon=mean(cluster_members$Lon),lat=mean(cluster_members$Lat)),
                             maptype="roadmap",color="bw") 
  }
  if(min(x_limits) < min(background_layer$data$lon) | max(x_limits) > max(background_layer$data$lon) |
     min(y_limits) < min(background_layer$data$lat) | max(y_limits) > max(background_layer$data$lat)){
    while(!exists("background_layer")){
      background_layer <- qmap(location = c(lon=mean(cluster_members$Lon),lat=mean(cluster_members$Lat)),
                               maptype="roadmap",color="bw",zoom=9) 
    }
  }
  
  gg <- background_layer
  gg <- gg + scale_x_continuous(limits = x_limits, expand = c(0, 0))
  gg <- gg + scale_y_continuous(limits = y_limits, expand = c(0, 0))
  gg <- gg + geom_polygon(data = myCircles_clu, aes(lon, lat, group = ID), color = NA, fill="red",alpha = 0.2)
  gg <- gg + geom_point(data=cluster_members, aes(x=Lon, y=Lat), color="black", shape=4, size=.5, stroke=0.5)
  gg <- gg + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color="yellow",size=0.2,data = coords)
  gg <- gg + coord_fixed(1.3)
  ggsave(paste0("maps/map_",d,"_",s,"_",i,".png",sep=""),gg)
  gg_net <- ggnet2(simplify(links),node.size = 3, node.color = "black", edge.size = 0.5, edge.color = "grey",mode = "fruchtermanreingold")
  ggsave(paste0("networks/net_",d,"_",s,"_",i,".png",sep=""),gg_net)
  rm(background_layer)
  
  progress(p,length(unique(inv_loc_clusters$Cluster_ID)))
  Sys.sleep(0.01)
  if (p == length(unique(inv_loc_clusters$Cluster_ID))) cat("\nClusters are mapped!")
  p <- p+1
}

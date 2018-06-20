# Packages needed for the cluster identification
# install.packages(c("data.table","igraph","prodlim","geosphere","svMisc"))

library(data.table)
library(igraph)
library(prodlim)
library(geosphere)
library(svMisc)

# Packages needed for the cluster visualization
# install.packages(c("ggplot2","ggmap","GGally))
library(ggplot2)
library(ggmap)
library(GGally)

options(stringsAsFactors = FALSE)

setwd("~")
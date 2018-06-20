source("global.R")

invpat_raw <- read.csv("data/invpat.csv")

## NB: the original disambiguation algorithm provides two strings with disambiguated codes
##     of inventor per each inventor-patent pair: "lower" and "upper".
##     It's up to researcher which one to use in the further analysis.

names(invpat_raw)[match("lower",names(invpat_raw))] <- "Invcode"

## drop unused columns of the data
invpat_raw <- invpat_raw[,c("Patent","Invcode","InvSeq","AppYearStr","Lat","Lon","Country","Class")]
gc()

## leave only the major class part of the class string
invpat_raw$Class <- sub("/.*","",invpat_raw$Class)
  
## drop duplicated inventor-patent pairs
invpat_raw <- invpat_raw[!duplicated(invpat_raw[,c("Invcode","Patent")]) & !duplicated(invpat_raw[,c("Invcode","Patent")],fromLast=T),]

## keep only patents with the first inventor from the US
invpat_raw <- data.table(invpat_raw)
invpat_raw <- invpat_raw[,FirstInvUS:=(Country[which.min(InvSeq)]=="US"),by=c("Patent")]
invpat_raw <- invpat_raw[invpat_raw$FirstInvUS,]

## drop inventor-patent pairs that are not geo-coded
invpat_raw <- invpat_raw[!is.na(invpat_raw$Lat) & !is.na(invpat_raw$Lon),]

## drop inventors with a non-unique location in a year of patent application
invpat_raw$LocCode <- as.numeric(factor(paste(invpat_raw$Lat,invpat_raw$Lon)))
invpat_raw <- invpat_raw[,UniqueLoc:=(length(unique(LocCode))==1),by=c("Invcode","AppYearStr")]
invpat_raw <- invpat_raw[invpat_raw$UniqueLoc,]

invpat_raw <- as.data.frame(invpat_raw)
save(invpat_raw,file="data/invpat.RData")
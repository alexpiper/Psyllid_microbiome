library(sp)
library(vegan)


# Bespoke function to convert a degree minute second coordinate into a fully numeric one,assuming D is used for degree sign

convert<-function(coord){
   t1 <- strsplit(coord, "D")
  d <- as.numeric(unlist(lapply(t1, "[", 1)))
  min <- as.numeric(unlist(sapply(strsplit(unlist(lapply(t1, "[",
                                                         2)),"'"), "[", 1)))
  sec <- as.numeric(substr(unlist(sapply(strsplit(unlist(lapply(t1,
                                                                "[", 2)),"'"), "[", 2)),1,5))
  return(d+min/60+sec/(60*60))
}

#Get the spatial distance between samples:
  envData <- read.csv("~/Google Drive/Martoni_Feb_2017_Analysis/Psyllids_data.csv", stringsAsFactors= FALSE)

  #Convert the locations into S and E

envData$S <- convert(unlist(lapply(strsplit(envData$Location," "),
                                   "[", 1)))
envData$E <- convert(unlist(lapply(strsplit(envData$Location," "),
                                   "[", 2)))

#Create a spatial distance between all samples
spatDist <- spDists(cbind(envData$E,envData$S)[!is.na(envData$E),],
                    longlat=TRUE)
rownames(spatDist) <- colnames(spatDist) <-
  envData$ID[!is.na(envData$E)]

#Get phylogenetic distance:
  phyloDist <- read.csv("~/Google
Drive/Martoni_Feb_2017_Analysis/Matrix1_phyloDist_Martoni.csv",
                        header=FALSE, stringsAsFactors=FALSE)
phyloSnum <- unlist(sapply(strsplit(phyloDist[,1], "_"),function(x)
  paste(x[1],x[2],sep="_")))

phyloDistM <- as.matrix(phyloDist[,-1],
                        dimnames=list(phyloSnum,phyloSnum))
colnames(phyloDistM) <- rownames(phyloDistM) <- phyloSnum

#Given two distance matrices, need to pair rows based on matching names.
intersectNames <- intersect(colnames(spatDist),colnames(phyloDistM))
phyloDistD <- as.dist(phyloDistM[match(intersectNames,
                                       rownames(phyloDistM)),match(intersectNames, colnames(phyloDistM))])
spatDistD <- as.dist(spatDist[match(intersectNames,
                                    rownames(spatDist)),match(intersectNames, colnames(spatDist))])
plot(phyloDistD~c(spatDistD+1),log="x")



#Add the molecular data from a previously saved file
load("~/Google Drive/Martoni_Feb_2017_Analysis/Psyllid_NGS_2017-03-08_small")
psyllidBacteriaOtus$length <- as.numeric(psyllidBacteriaOtus$length
)
psyllidBacteriaOtus$identity <-
  as.numeric(psyllidBacteriaOtus$identity )
rownames(communityM) <- gsub("-","_",rownames(communityM))
communityM <- communityM[rownames(communityM) %in% phyloSnum,]
plot(rowSums(communityM),log="y",main="Sequencing depth variability")

#NMDS
mds <- metaMDS(communityM, trymax=100)
#fails to converge with 100 tries.
Try with 3 axes:
  mds <- metaMDS(communityM, trymax=100, k=3)

#still fails to converge.
#Try with only samples > 500 sequence and only OTUs found in > 1 psyllid
mds <- metaMDS(communityM[rowSums(communityM)>500,colSums(communityM>0)>1],
          trymax=100)
plot(mds$points, cex=0)
text(mds$points,rownames(communityM[rowSums(communityM)>500,]),cex=0.5)
quartz(height=6,width=11) ##Use X11 on PC?
par(mfrow=c(1,2))

#Rank abundance graph:
  plot(sort(colSums(communityM),decreasing=TRUE),log="y")

#Frequency
plot(sort(colSums(communityM>0),decreasing=TRUE),log="y")

#what are most frequent:

  head(sort(colSums(communityM>0),decreasing=TRUE),20)
  
#get names: (This just takes the line above and uses match to find it in the OTU data)
psyllidBacteriaOtus[match(names(head(sort(colSums(communityM>0), decreasing=TRUE),20)), psyllidBacteriaOtus$otu),]
communityD <-as.matrix(vegdist(communityM[rowSums(communityM)>500, colSums(communityM>0)>1]))
intersectNames <- intersect(colnames(communityD),colnames(phyloDistM))

phyloDistD <- as.dist(phyloDistM[match(intersectNames, rownames(phyloDistM)),match(intersectNames, colnames(phyloDistM))])
communityDD <- as.dist(communityD[match(intersectNames, rownames(communityD)),match(intersectNames, colnames(communityD))])

plot(communityDD~phyloDistD)
mantel(communityDD, phyloDistD)
intersectNames <- intersect(colnames(communityD),colnames(spatDist))
spatDistD <- as.dist(spatDist[match(intersectNames,
                                    rownames(spatDist)),match(intersectNames, colnames(spatDist))])
communityDD <- as.dist(communityD[match(intersectNames,
                                        rownames(communityD)),match(intersectNames, colnames(communityD))])
plot(communityDD~c(spatDistD+1), log="x")
mantel(communityDD, spatDistD)
Three way intersect and mantel tests
intersectNames <-
  intersect(intersect(colnames(communityD),colnames(spatDist)),colname
            s(phyloDistM))
spatDistD <- as.dist(spatDist[match(intersectNames,
                                    rownames(spatDist)),match(intersectNames, colnames(spatDist))])
communityDD <- as.dist(communityD[match(intersectNames,
                                        rownames(communityD)),match(intersectNames, colnames(communityD))])
phyloDistD <- as.dist(phyloDistM[match(intersectNames,
                                       rownames(phyloDistM)),match(intersectNames, colnames(phyloDistM))])
mantel.partial(communityDD, spatDistD, phyloDistD)
mantel.partial(communityDD, log(spatDistD+1), phyloDistD) ##Log
distance not better
mantel.partial(communityDD, phyloDistD, spatDistD)

#Adonis approach 

#Only look at species with > 2 records
speciesFreq <- table(envData$Species[match(rownames(communityD),envData$ID)])
subCommunity <- communityM[rowSums(communityM)>500,colSums(communityM>0)>1]
freqSppSamples <- envData$ID[envData$Species %in%
                               names(speciesFreq)[speciesFreq>2]]
subCommunity <- subCommunity[rownames(subCommunity) %in%
                               freqSppSamples,]
speciesSubCom <- envData$Species[match(rownames(subCommunity),
                                       envData$ID)]
adonis(subCommunity~as.factor(speciesSubCom))

#Clustering:
  plot(hclust(vegdist(communityM[rowSums(communityM)>500,colSums(communityM>0)>0])))

#Add plant phylogenetic distance
test <- read.table("/Volumes/NO NAME/00PLANTS/NGSplantcorrect.txt",
                   header=TRUE, sep="\t")
test <- test[,-1]
rownames(test) <- colnames(test)
speciesList <- unlist(sapply(strsplit(rownames(test), "\\."), "[",
                             1))
for(i in ncol(test))
{
  test[speciesList == speciesList[i],i] <- 0
  
  
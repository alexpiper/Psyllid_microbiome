run_mantel <- function(matrix, dists, samples, type="mantel"){
  if(type=="mantel" && class(dists)=="character"){
    dists <- enframe(dists) %>% dplyr::rename(dist1 = value)
  } else if (type=="partial" && class(dists)=="character"){
    message("Expanding all possible combinations for partial mantel")
    dists <- expand_grid(dist1 = dists, dist2=dists) %>%
      filter(!dist1==dist2)
  }
  if(type=="mantel"){
    out <- dists %>% 
      split(rownames(.)) %>%
      purrr::map(function(y){
        as.data.frame(mantel(matrix[samples, samples],
                             get(y$dist1)[samples, samples]
        )[c("statistic","signif", "permutations")]) %>%
          mutate(dist1 = y$dist1, type="mantel") 
      })%>%
      bind_rows()
  }else if (type=="partial"){
    out <- dists %>% 
      split(rownames(.)) %>%
      purrr::map(function(y){
        as.data.frame(mantel.partial(matrix[samples, samples],
                                     get(y$dist1)[samples, samples], 
                                     get(y$dist2)[samples,samples]
        )[c("statistic","signif", "permutations")]) %>%
          mutate(dist1 = y$dist1, dist2=y$dist2, type="partial_mantel") 
      }) %>%
      bind_rows()
  }
  return(out)
}

# RFmeasures --------------------------------------------------------------
tree_measures <- function(tree1, tree2, measure="all", runs=99){
  if(!measure %in% c("all","RobinsonFoulds","PhylogeneticInfoDistance",
                     "ClusteringInfoDist", "NyeTreeDist",
                     "MatchingSplitInfoDistance")){
    stop("Error, invalid distance")
  }
    
  if (measure=="all"){
    measure <- c("RobinsonFoulds","PhylogeneticInfoDistance",
                 "ClusteringInfoDist", "NyeTreeDist",
                 "MatchingSplitInfoDistance")
  }
  #transform trees
  tree1 <- as.phylo(tree1)
  tree2 <- as.phylo(tree2)
  if(is.rooted(tree1)){
    cat("tree1 is rooted. Unrooting tree1.\n")
    tree1 <- unroot(tree1)
  }
  if(is.rooted(tree2)){
    cat("tree2 is rooted. Unrooting tree2.\n")
    tree2 <- unroot(tree2)
  }
  
  # Create null models
  #Shuffle tips ("null model")
  nulltree <- vector("list", length=1)
  nulltree[[1]] <- tree2
  nulltree <- lapply(rep(nulltree, runs), tipShuffle)
  class(nulltree) <- "multiPhylo"
  
  #new (random) trees ("random model")
  tips <- as.phylo(tree2)$tip.label
  randomtree <- rmtree(N=runs,n=length(tips),tip.label = tips)

  #Create lists to store objects
  tmp <- vector("list", length=length(measure))
  names(tmp) <- measure
  #Process observed trees
  if("RobinsonFoulds" %in% measure){
    obs  <- RobinsonFoulds(tree1, tree2, normalize=TRUE)
    null <- sapply(1:runs, FUN=function(x){RobinsonFoulds(tree1, nulltree[[x]], normalize=TRUE)})
    rd <- sapply(1:runs, FUN=function(x){RobinsonFoulds(tree1, randomtree[[x]], normalize=TRUE)})
    null.mean <- mean(null)
    null.sd <- sd(null)
    obs.rank.null <- apply(X = rbind(obs, null), MARGIN = 2, FUN = rank)[1, ]
    obs.rank.null <- ifelse(is.na(null.mean), NA, obs.rank)
    rd.mean <- mean(rd)
    rd.sd <- sd(rd)
    obs.rank.rd <- apply(X = rbind(obs,rd), MARGIN = 2, FUN = rank)[1, ]
    obs.rank.rd <- ifelse(is.na(rd.mean), NA, obs.rank)
    
    tmp$RobinsonFoulds <- data.frame(obs=obs, obs.rank.null = obs.rank.null,
                                     null.mean = null.mean, null.sd = null.sd ,
                                     rd.mean = rd.mean, rd.sd = rd.sd,
                                     obs.z = 
                                     obs.p = obs.rank/(runs + ),
                                     pvalRd = sum(obs > rd)/(runs + 1)) 
    
    
  }
  if("PhylogeneticInfoDistance" %in% measure){
    obs  <- PhylogeneticInfoDistance(tree1, tree2, normalize=TRUE)
    null <- sapply(1:runs, FUN=function(x){PhylogeneticInfoDistance(tree1, nulltree[[x]], normalize=TRUE)})
    rd <- sapply(1:runs, FUN=function(x){PhylogeneticInfoDistance(tree1, randomtree[[x]], normalize=TRUE)})
    tmp$PhylogeneticInfoDistance <- data.frame(stat=obs, pval=sum(obs > null)/(runs + 1), pvalRd = sum(obs > rd)/(runs + 1)) 
  }
  if("ClusteringInfoDist" %in% measure){
    obs  <- ClusteringInfoDist(tree1, tree2, normalize=TRUE)
    null <- sapply(1:runs, FUN=function(x){ClusteringInfoDist(tree1, nulltree[[x]], normalize=TRUE)})
    rd <- sapply(1:runs, FUN=function(x){ClusteringInfoDist(tree1, randomtree[[x]], normalize=TRUE)})
    tmp$ClusteringInfoDist <- data.frame(stat=obs, pval=sum(obs > null)/(runs + 1), pvalRd = sum(obs > rd)/(runs + 1)) 
  }
  if("NyeTreeDist" %in% measure){
    obs  <- 1-NyeTreeSimilarity(tree1, tree2, normalize=TRUE)
    null <- 1-sapply(1:runs, FUN=function(x){NyeTreeSimilarity(tree1, nulltree[[x]], normalize=TRUE)})
    rd <- 1-sapply(1:runs, FUN=function(x){NyeTreeSimilarity(tree1, randomtree[[x]], normalize=TRUE)})
    tmp$NyeTreeDist <- data.frame(stat=obs, pval=sum(obs > null)/(runs + 1), pvalRd = sum(obs > rd)/(runs + 1)) 
  }
  if("MatchingSplitInfoDistance" %in% measure){
    obs  <- MatchingSplitInfoDistance(tree1, tree2, normalize=TRUE)
    null <- sapply(1:runs, FUN=function(x){MatchingSplitInfoDistance(tree1, nulltree[[x]], normalize=TRUE)})
    rd <- sapply(1:runs, FUN=function(x){MatchingSplitInfoDistance(tree1, randomtree[[x]], normalize=TRUE)})
    tmp$MatchingSplitInfoDistance <- data.frame(stat=obs, pval=sum(obs > null)/(runs + 1), pvalRd = sum(obs > rd)/(runs + 1)) 
  }
  out <- bind_rows(tmp, .id="measure")
  return(out)
}


## Try treedist metrics
tree1  <- hclust(as.dist(spat.dist[samples,samples]), method = "ward.D2") # Average = average
tree2 <- hclust(as.dist(distlist$Aitchison[samples,samples]), method = "ward.D2")
runs=100
set.seed(999)

tree1 <- as.phylo(tree1)
tree2 <- as.phylo(tree2)
VisualizeMatching(MatchingSplitInfoDistance,tree1, tree2)

tree_measures(tree1, tree2, runs=99)

# Test with a random tree

nulltree <- shuffleTips(1, tree2)
tree_measures(tree1, nulltree, runs=99)
library('TreeDist')




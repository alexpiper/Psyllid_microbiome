
# Run mantel vegan --------------------------------------------------------------

run_mantel_vegan <- function(x, dists, subsample, type="mantel"){
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
        as.data.frame(vegan::mantel(x[subsample, subsample],
                             get(y$dist1)[subsample, subsample]
        )[c("statistic","signif", "permutations")]) %>%
          mutate(dist1 = y$dist1, type="mantel") 
      })%>%
      bind_rows()
  }else if (type=="partial"){
    out <- dists %>% 
      split(rownames(.)) %>%
      purrr::map(function(y){
        as.data.frame(vegan::mantel.partial(x[subsample, subsample],
                                     get(y$dist1)[subsample, subsample], 
                                     get(y$dist2)[subsample,subsample]
        )[c("statistic","signif", "permutations")]) %>%
          mutate(dist1 = y$dist1, dist2=y$dist2, type="partial_mantel") 
      }) %>%
      bind_rows()
  }
  return(out)
}

# Permustats
#test <- vegan::mantel(veg.dist, env.dist)
#perm <- vegan::permustats(test)
#test3 <- summary(perm)

# Run mantel ecodist ------------------------------------------------------
# need to generalize to take in N matrices rather than just 3 
run_mantel <- function(x, dists, subsample, type="mantel", nboot=1000){
  if(type=="mantel" && class(dists)=="character"){
    dists <- enframe(dists) %>% dplyr::rename(dist1 = value)
  } else if (type=="partial" && class(dists)=="character"){
    message("Expanding all possible combinations for partial mantel")
    dists <- combinat::permn(dists) %>%
      purrr::map(t) %>%
      purrr::map(as_tibble) %>%
      dplyr::bind_rows() %>%
      dplyr::filter(!duplicated(V1))
  }
  if(type=="mantel"){
    out <- dists %>% 
      split(rownames(.)) %>%
      purrr::map(function(y){
        x_mat <- x[subsample, subsample]
        y_mat <- get(y$dist1)[subsample, subsample]
        as.data.frame(t(ecodist::mantel( lower(x_mat) ~ lower(y_mat), nboot=nboot))) %>%
          mutate(dist1 = y$dist1, type="mantel", nboot=1000) 
      })%>%
      bind_rows()
  }else if (type=="partial"){
    out <- dists %>% 
      split(rownames(.)) %>%
      purrr::map(function(y){
        x_mat <- x[subsample, subsample]
        y_mat <- get(y$V1)[subsample, subsample]
        z1_mat <- get(y$V2)[subsample, subsample]
        z2_mat <- get(y$V3)[subsample, subsample]
        
        as.data.frame(t(ecodist::mantel(lower(x_mat) ~ lower(y_mat) + lower(z1_mat) + lower(z2_mat), nboot=nboot))) %>%
          mutate(dist1 = y$V1, dist2=y$V2, dist3=y$V3, type="partial_mantel") 
      }) %>%
      bind_rows()
  }
  return(out)
}

# compare_trees --------------------------------------------------------------
compare_trees <- function(tree1, tree2, measure="all", runs=99){
  availmeasures <- c("RobinsonFoulds","PhylogeneticInfoDistance",
                    "ClusteringInfoDist", "NyeTreeDist",
                    "MatchingSplitInfoDistance")
  if (length (measure)== 1 && measure=="all"){
    measure <- availmeasures
  }
  if(any(!measure %in% availmeasures)){
    stop("Error, invalid TreeDist measure :, ", measure[!measure %in% availmeasures])
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
  # Create null model from tree2 by shuffling tips ("taxa.labels")
  nulltree <- vector("list", length=1)
  nulltree[[1]] <- tree2
  nulltree <- lapply(rep(nulltree, runs), tipShuffle)
  class(nulltree) <- "multiPhylo"

  #Create lists to store objects
  tmp <- vector("list", length=length(measure))
  names(tmp) <- measure

  #Calculate metrics
  if("RobinsonFoulds" %in% measure){
    obs  <- RobinsonFoulds(tree1, tree2, normalize=TRUE)
    null <- sapply(1:runs, FUN=function(x){RobinsonFoulds(tree1, nulltree[[x]], normalize=TRUE)})
    null.mean <- mean(null)
    null.sd <- sd(null)
    obs.rank <- apply(X = rbind(obs, null), MARGIN = 2, FUN = rank)[1, ]
    obs.rank <- ifelse(is.na(null.mean), NA, obs.rank)
    tmp$RobinsonFoulds <- data.frame(obs=obs, obs.rank = obs.rank,
                                     null.mean = null.mean, null.sd = null.sd,
                                     obs.z = (obs - null.mean)/null.sd,
                                     obs.p = obs.rank/(runs + 1 )) 
    }
  if("PhylogeneticInfoDistance" %in% measure){
    obs  <- PhylogeneticInfoDistance(tree1, tree2, normalize=TRUE)
    null <- sapply(1:runs, FUN=function(x){PhylogeneticInfoDistance(tree1, nulltree[[x]], normalize=TRUE)})
    null.mean <- mean(null)
    null.sd <- sd(null)
    obs.rank <- apply(X = rbind(obs, null), MARGIN = 2, FUN = rank)[1, ]
    obs.rank <- ifelse(is.na(null.mean), NA, obs.rank)
    tmp$PhylogeneticInfoDistance <- data.frame(obs=obs, obs.rank = obs.rank,
                                     null.mean = null.mean, null.sd = null.sd,
                                     obs.z = (obs - null.mean)/null.sd,
                                     obs.p = obs.rank/(runs + 1 )) 
  }
  if("ClusteringInfoDist" %in% measure){
    obs  <- ClusteringInfoDist(tree1, tree2, normalize=TRUE)
    null <- sapply(1:runs, FUN=function(x){ClusteringInfoDist(tree1, nulltree[[x]], normalize=TRUE)})
    null.mean <- mean(null)
    null.sd <- sd(null)
    obs.rank <- apply(X = rbind(obs, null), MARGIN = 2, FUN = rank)[1, ]
    obs.rank <- ifelse(is.na(null.mean), NA, obs.rank)
    tmp$ClusteringInfoDist <- data.frame(obs=obs, obs.rank = obs.rank,
                                               null.mean = null.mean, null.sd = null.sd,
                                               obs.z = (obs - null.mean)/null.sd,
                                               obs.p = obs.rank/(runs + 1 )) 
  }
  if("NyeTreeDist" %in% measure){
    obs  <- 1-NyeTreeSimilarity(tree1, tree2, normalize=TRUE)
    null <- 1-sapply(1:runs, FUN=function(x){NyeTreeSimilarity(tree1, nulltree[[x]], normalize=TRUE)})
    null.mean <- mean(null)
    null.sd <- sd(null)
    obs.rank <- apply(X = rbind(obs, null), MARGIN = 2, FUN = rank)[1, ]
    obs.rank <- ifelse(is.na(null.mean), NA, obs.rank)
    tmp$NyeTreeDist <- data.frame(obs=obs, obs.rank = obs.rank,
                                         null.mean = null.mean, null.sd = null.sd,
                                         obs.z = (obs - null.mean)/null.sd,
                                         obs.p = obs.rank/(runs + 1 )) 
  }
  if("MatchingSplitInfoDistance" %in% measure){
    obs  <- MatchingSplitInfoDistance(tree1, tree2, normalize=TRUE)
    null <- sapply(1:runs, FUN=function(x){MatchingSplitInfoDistance(tree1, nulltree[[x]], normalize=TRUE)})
    null.mean <- mean(null)
    null.sd <- sd(null)
    obs.rank <- apply(X = rbind(obs, null), MARGIN = 2, FUN = rank)[1, ]
    obs.rank <- ifelse(is.na(null.mean), NA, obs.rank)
    tmp$MatchingSplitInfoDistance <- data.frame(obs=obs, obs.rank = obs.rank,
                                  null.mean = null.mean, null.sd = null.sd,
                                  obs.z = (obs - null.mean)/null.sd,
                                  obs.p = obs.rank/(runs + 1 )) 
  }
  out <- bind_rows(tmp, .id="measure")
  return(out)
}

# Run tree comparison  --------------------------------------------------------------
run_tree_comparison <- function(x, dists, subsample, clust.method = "average", measure="all", runs = 100){
    dists <- enframe(dists) %>% dplyr::rename(dist1 = value)
    out <- dists %>% 
      split(rownames(.)) %>%
      purrr::map(function(y){
        tree1  <- as.phylo(hclust(as.dist(get(y$dist1)[subsample, subsample]), method = clust.method))
        tree2 <- as.phylo(hclust(as.dist(x[subsample, subsample]), method = clust.method))
        compare_trees(tree1, tree2, measure=measure, runs=runs) %>%
          dplyr::mutate(dist = y$dist1)
      })%>%
      bind_rows()
  return(out)
}


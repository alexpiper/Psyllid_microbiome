# Negative to positive eig ------------------------------------------------
neg_to_pos_eig = function(PCo){
  min_val = PCo %>% as.vector %>% min
  print(min_val)
  if(min_val < 0){
    PCo = PCo + abs(min_val)
  }
  return(PCo)
}


# Average descendents -----------------------------------------------------
# average a ggtree treedata value from tips over all nodes
average_descendants <- function(x, tree, df) {
  descendants <- getDescendants(tree, x)
  df %>%
    dplyr::filter(node %in% descendants) %>%
    dplyr::mutate(values = case_when(
      !is.na(values) ~ values,
      is.na(values) ~ 0
    )) %>%
    pull(values) %>%
    mean(na.rm=TRUE)
}


# Pairwise adonis ---------------------------------------------------------

# From https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis2.R
pairwise.adonis2 <- function(x, data, strata = NULL, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis(xnew,data=mdat1, ... )
    }else{ad <- adonis(xnew,data=mdat1,strata= mdat1[,strata], ... )}
    
    res[nameres[elem+1]] <- ad[1]
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 


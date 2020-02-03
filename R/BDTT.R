

# New to old OTUs ---------------------------------------------------------

New_toOld_OTUs_castor <- function(similarity,tree){
  Ntips=length(tree$tip.label)
  NameTips=tree$tip.label
  
  #Collapse tips accorindg to the resolution 
  ColapsedNodes=collapse_tree_at_resolution(tree,similarity)
  DescendingTips=lapply(ColapsedNodes$collapsed_nodes,tree,FUN=function(i,tree){get_subtree_at_node(tree, i)$subtree$tip.label})
  names(DescendingTips)=as.character(ColapsedNodes$collapsed_nodes+Ntips)
  
  #Lenght of the different categories of OTUs
  N_newOTUS=length(DescendingTips)
  N_collapsedOldOtus= sum(do.call(rbind, lapply(DescendingTips, function(x) length(x))))
  N_totalNewOTUS=Ntips-N_collapsedOldOtus+N_newOTUS
  
  print(paste(similarity," similarity provides ",N_totalNewOTUS," total new OTUs",sep=""))
  
  NewOTUs_OldOTUs_Matrix=NA
  if (N_totalNewOTUS>2) {
    #Names of the different categories of OTUs
    collapsedOldOtus= unlist(DescendingTips)
    UncollapsedOTUs=NameTips[!NameTips%in%collapsedOldOtus]
    
    #Creating New to Old matrix correpondance
    NewOTUSNames=c(UncollapsedOTUs,names(DescendingTips))
    NewOTUs_OldOTUs_Matrix=Matrix(0, ncol=N_totalNewOTUS, nrow=Ntips, dimnames = list(tree$tip.label,NewOTUSNames))
    
    #Filling  New to Old matrix correpondance / uncollapsed tips
    if (length(UncollapsedOTUs)>1){diag(NewOTUs_OldOTUs_Matrix[UncollapsedOTUs,UncollapsedOTUs])=1}     #just fill the diagonal as tips are unchanged
    if (length(UncollapsedOTUs)==1){NewOTUs_OldOTUs_Matrix[UncollapsedOTUs,UncollapsedOTUs]=1}
    
    #Filling  New to Old matrix correpondance / collapsed tips
    for (i in names(DescendingTips)){NewOTUs_OldOTUs_Matrix[DescendingTips[[i]],i]=1}
    
  } else{print("Only 2 OTUS (or less) present at this resolution -- are you sure this is meaningfull?")} 
  
  newtree <- ColapsedNodes$tree
  newtree$tip.label <- NewOTUSNames
  out=list(CorresMatrix = NewOTUs_OldOTUs_Matrix, tree=newtree)
  return(out)
}

getNew_OTUs_sample_matrix=function(similarity,sampleOTUs,tree){
  New2old <- New_toOld_OTUs_castor(similarity=similarity,tree=tree) 
  CorresMatrix <- New2old$CorresMatrix
  newtree <- New2old$tree
  
  sample=Matrix(t(sampleOTUs))
  sp=colnames(sample)
  Newsample=t(sample[,sp] %*% CorresMatrix[sp,]) # multiplies them
  out=list(NewSample=as.matrix(Newsample),CorrepondanceOTUs=CorresMatrix, tree=newtree)
  return(out)
}


# Get Beta diversity ------------------------------------------------------

getBDTT <- function(similarity, tree, sampleOTUs, metrics=c("Bray","Jaccard","Aitchison","Philr"), zeroes="Impute", parallel=FALSE, onlyBeta=TRUE, quiet=FALSE){
  New_OTUs_sample_matrix=getNew_OTUs_sample_matrix(similarity=similarity, sampleOTUs=sampleOTUs, tree=tree)
  OTU_table=t(New_OTUs_sample_matrix$NewSample)
  newtree <- New_OTUs_sample_matrix$tree
  
  # If compositional stats - Deal with zeroes
  if (any(c("Aitchison","Philr") %in% metrics) & stringr::str_to_lower(zeroes) =="pseudo"){
    if(!quiet) {message("Adding a pseudocount for compositional distances")}
    OTU_table_n0 <- prop.table(OTU_table + 1, margin=1)
  } else if (any(c("Aitchison","Philr") %in% metrics) & stringr::str_to_lower(zeroes) =="impute"){
    if(!quiet) {message("Imputing zeroes for compositional distances")}
    OTU_table_n0 <- as.matrix(zCompositions::cmultRepl(OTU_table, method="CZM", output="prop", suppress.print=quiet))
  }
  
  # If phylogenetic distances are used - root tree
  if (any(c("Philr","Unifrac","WUnifrac") %in% metrics)){
    newtree  <- multi2di(newtree)
  }
  ## Calculate distances
  # Bray curtis distance
  if ("Bray" %in% metrics){
    Bray <- as.matrix(vegdist(OTU_table, method="bray"))
  }
  
  #Jaccard distance
  if ("Jaccard" %in% metrics){
    Jaccard <- as.matrix(vegdist(OTU_table, method="jac",binary = T))
  }
  #JacTT=as.matrix(designdist(OTU_table, "(2*pmin(b,c))/(a+2*pmin(b,c))",abcd=TRUE,terms = "binary"))
  
  # Aitchison distance
  if ("Aitchison" %in% metrics){
    Aitchison <- as.matrix(vegdist(CoDaSeq::codaSeq.clr(OTU_table_n0), method="euclidean")) # Replace this with custom CLR function
  }
  
  #PhilR distance 
  if ("Philr" %in% metrics){
    if(!quiet){
      message("Calculating PhilR distance")
    Philr <- as.matrix(vegdist(philr::philr(OTU_table_n0, newtree,
                       part.weights='enorm.x.gm.counts',
                       ilr.weights='blw.sqrt'), method="euclidean"))
    message("PhilR distance complete")
    } else if(quiet){
      suppressMessages(
        Philr <- as.matrix(vegdist(philr::philr(OTU_table_n0, newtree,
                                              part.weights='enorm.x.gm.counts',
                                              ilr.weights='blw.sqrt'), method="euclidean"))
      )
    }
  }
  
  # Unweighted Unifrac distance
  if ("Unifrac" %in% metrics){
    if(!quiet) {message("Calculating Unweighted Unifrac distance...")}
    physeq <- phyloseq(otu_table(OTU_table, taxa_are_rows=FALSE), phy_tree(newtree))
    Unifrac <- as.matrix(phyloseq::UniFrac(physeq, weighted=FALSE, parallel = parallel))
    if(!quiet) {message("Unweighted Unifrac distance complete")}
  }
  
  # Weighted Unifrac distance
  if ("WUnifrac" %in% metrics){
    if(!quiet) {message("Calculating Weighted Unifrac distance...")}
    physeq <- phyloseq(otu_table(OTU_table, taxa_are_rows=FALSE), phy_tree(newtree))
    WUnifrac <- as.matrix(phyloseq::UniFrac(physeq, weighted=TRUE, parallel = parallel))
    if(!quiet) {message("Unweighted Unifrac distance complete")}
  }
  
  # Abind only those that are in metrics
  AllBetas=abind(lapply(metrics, get, envir=sys.frame(sys.parent(0))), along = 0)
  dimnames(AllBetas)[[1]]=metrics
  
  if (onlyBeta==FALSE){res <- list(Beta_Div=AllBetas, 
                            NewOTU_table=New_OTUs_sample_matrix$NewSample, 
                            New_to_old_OTUs=New_OTUs_sample_matrix$CorrepondanceOTU,
                            Newtree = New_OTUs_sample_matrix$tree )}
  if (onlyBeta==TRUE){res <- AllBetas}
  
  return(res)
}


# Beta diversity function -------------------------------------------------


BDTT <- function(similarity_slices, tree, sampleOTUs, metrics=c("Bray","Jaccard","Aitchison","Philr"), zeroes="Impute", parallel=FALSE, onlyBeta=TRUE, quiet=FALSE){
  
  #Check orientation of samples
  if(!taxa_are_rows(sampleOTUs)){ 
    sampleOTUs <- t(sampleOTUs)
    }
  
  Betas <- base::lapply(similarity_slices, getBDTT, tree=tree, sampleOTUs=sampleOTUs, metrics=metrics, zeroes=zeroes, parallel=parallel, onlyBeta=onlyBeta, quiet=quiet)
  names(Betas) <- similarity_slices
  out <- do.call(function(...){abind(...,along=0)}, Betas)
  return(out)
}



# Sloans neutral model through time ---------------------------------------


NMTT <- function(similarity_slices, tree, sampleOTUs, return="summary", quiet=FALSE){
  #Check orientation of samples
  if(!taxa_are_rows(sampleOTUs)){ 
    sampleOTUs <- t(sampleOTUs)
  }
  
  getNMTT <- function(similarity, tree, sampleOTUs, quiet=FALSE){
    New_OTUs_sample_matrix <- getNew_OTUs_sample_matrix(similarity=similarity, sampleOTUs=sampleOTUs, tree=tree)
    OTU_table <- t(New_OTUs_sample_matrix$NewSample)
    newtree <- New_OTUs_sample_matrix$tree
    out <- fit_sncm(OTU_table, pool=OTU_table)
    #out$tree <- newtree
    #out$otu <- OTU_table
    return(out)
  }
  
  neutral <- lapply(similarity_slices, getNMTT, tree=tree, sampleOTUs=sampleOTUs, quiet=quiet)
  names(neutral) <- similarity_slices
  
  # Create plots
  plotdata <- neutral %>%
    purrr::map_df("predictions", .id="slice") %>%
    left_join(purrr::map_df(neutral, "fitstats", .id="slice"), by="slice")
  
  plots <- ggplot(plotdata) + 
    geom_point(aes(x = log(p), y = freq, fill = fit_class),
               shape = 21, color = "black", size = 2, alpha = 0.75) +
    geom_text(aes(x=mean(log(p)), label = paste("r^2 ==", round(Rsqr, 4)), y = 0.95, size = 5 ))+
    geom_text(aes(x=mean(log(p)), label = paste("m ==", round(m,  4)), y = 0.9, size = 5 ))+
    geom_line(aes(x = log(p), y = freq.pred), color = "blue") + 
    geom_line(aes(x = log(p), y = pred.lwr), color = "blue", linetype = "dashed") + 
    geom_line(aes(x = log(p), y = pred.upr), color = "blue", linetype = "dashed") +
    theme_bw() +
    facet_wrap(~slice)+
    xlab("log(Mean Relative Abundance)") +
    ylab("Frequency")
  
  message("made plots")
  if(return == "summary"){
    out <- list()
    out$summary <- neutral %>%
      purrr::map_df("predictions", .id="slice") %>%
      group_by(slice, fit_class)%>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      ungroup() %>%
      mutate(fit_class = fit_class %>%
               str_replace("\\ (.*)$", "") %>%
               str_replace("As", "Predicted")) %>%
      pivot_wider(id_cols=slice, 
                  names_from = fit_class, 
                  values_from = c(n, freq)) %>%
      left_join(purrr::map_df(neutral, "fitstats", .id="slice"), by="slice")
    
    out$plots <- plots
  } else if (return == "all"){
    out <- neutral
  } else if (return == "plots"){
    out <- plots
  }
  return(out)
}

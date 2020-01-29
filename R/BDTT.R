
New_toOld_OTUs_castor=function(similarity,tree)
{
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
  if (N_totalNewOTUS>2) 
  {
    #Names of the different categories of OTUs
    collapsedOldOtus= unlist(DescendingTips)
    UncollapsedOTUs=NameTips[!NameTips%in%collapsedOldOtus]
    
    #Creating New to Old matrix correpondance
    NewOTUSNames=c(UncollapsedOTUs,names(DescendingTips))
    NewOTUs_OldOTUs_Matrix=Matrix(0,ncol=N_totalNewOTUS,nrow=Ntips,dimnames = list(tree$tip.label,NewOTUSNames))
    
    #Filling  New to Old matrix correpondance / uncollapsed tips
    if (length(UncollapsedOTUs)>1){diag(NewOTUs_OldOTUs_Matrix[UncollapsedOTUs,UncollapsedOTUs])=1}     #just fill the diagonal as tips are unchanged
    if (length(UncollapsedOTUs)==1){NewOTUs_OldOTUs_Matrix[UncollapsedOTUs,UncollapsedOTUs]=1}
    
    #Filling  New to Old matrix correpondance / collapsed tips
    for (i in names(DescendingTips)){NewOTUs_OldOTUs_Matrix[DescendingTips[[i]],i]=1}
    
  } else{print("Only 2 OTUS (or less) present at this resolution -- are you sure this is meaningfull?")} 
  
  
  return(NewOTUs_OldOTUs_Matrix)
}

getNew_OTUs_sample_matrix=function(similarity,sampleOTUs,tree)
{
  CorresMatrix=New_toOld_OTUs_castor(similarity=similarity,tree=tree) 
  sample=Matrix(t(sampleOTUs))
  sp=colnames(sample)
  Newsample=t(sample[,sp] %*% CorresMatrix[sp,]) # multiplies them
  out=list(NewSample=as.matrix(Newsample),CorrepondanceOTUs=CorresMatrix)
  return(out)
}

getBeta=function(mat,ab=F)
{
  if (ab==T)
  {   
    h=bray.part(data.frame(mat))
    bctu=as.matrix(h[[1]])
    bcne=as.matrix(h[[2]])
    bc=as.matrix(h[[3]])
    res=abind(bctu,bcne,bc,along=0)
    dimnames(res)[[1]]=c("bctu","bcne","bc")
  }
  if (ab==F)
  {
    mat[mat>0]=1
    h=beta.pair(data.frame(mat), index.family="jaccard")
    hh=beta.pair(data.frame(mat), index.family="sorensen")
    jtu=as.matrix(h[["beta.jtu"]])
    jne=as.matrix(h[["beta.jne"]])
    jac=as.matrix(h[["beta.jac"]])
    stu=as.matrix(hh[["beta.sim"]])
    sne=as.matrix(hh[["beta.sne"]])
    sor=as.matrix(hh[["beta.sor"]])  
    res=abind(jtu,jne,jac,stu,sne,sor,along=0)
    dimnames(res)[[1]]=c("jtu","jne","jac","stu","sne","sor")  
  }
  return(res)
}

#getBDTT=function(similarity,tree,sampleOTUs,onlyBeta=T,metric="jac")
#{
#  New_OTUs_sample_matrix=getNew_OTUs_sample_matrix(similarity=similarity,sampleOTUs=sampleOTUs,tree=tree)


#  BetasPA=getBeta(t(New_OTUs_sample_matrix[[1]]),ab=F)[c("jtu","jac"),,]
#rownames(BetasPA)=paste(rownames(BetasPA),"_Simi_",similarity,sep="")
#  BetasAb=getBeta(t(New_OTUs_sample_matrix[[1]]),ab=T)[c("bctu","bc"),,]  
#rownames(BetasAb)=paste(rownames(BetasAb),"_Simi_",similarity,sep="")
#  AllBetas=abind(BetasPA,BetasAb,along = 1)

#  AllBetas=AllBetas[metric,,]

#  if (onlyBeta==F){res=list(Beta_Div=AllBetas,NewOTU_table=New_OTUs_sample_matrix[[1]],New_to_old_OTUs=New_OTUs_sample_matrix[[2]])}
#  if (onlyBeta==T){res=AllBetas}

#  return(res)
#}


getBDTT=function(similarity,tree,sampleOTUs,onlyBeta=T)
{
  New_OTUs_sample_matrix=getNew_OTUs_sample_matrix(similarity=similarity,sampleOTUs=sampleOTUs,tree=tree)
  OTU_table=t(New_OTUs_sample_matrix[[1]])
  
  Bray=as.matrix(vegdist(OTU_table,method="bray"))
  Jac=as.matrix(vegdist(OTU_table,method="jac",binary = T))
  JacTT=as.matrix(designdist(OTU_table, "(2*pmin(b,c))/(a+2*pmin(b,c))",abcd=TRUE,terms = "binary"))
  
  #Compositional
  aitch_pseudo <- as.matrix(vegdist(codaSeq.clr(prop.table(OTU_table + 1, margin=1)), method="euclidean"))
  aitch_rep <- as.matrix(vegdist(codaSeq.clr(cmultRepl(OTU_table, method="CZM", output="prop")), method="euclidean"))
  
  
  AllBetas=abind(Bray,Jac,JacTT,aitch_pseudo,aitch_rep,along = 0)
  dimnames(AllBetas)[[1]]=c("Bray","Jac","Jac_TT","Aitch_pseudo","Aitch_rep")
  
  if (onlyBeta==F){res=list(Beta_Div=AllBetas,NewOTU_table=New_OTUs_sample_matrix[[1]],New_to_old_OTUs=New_OTUs_sample_matrix[[2]])}
  if (onlyBeta==T){res=AllBetas}
  
  return(res)
}

BDTT=function(similarity_slices,tree,sampleOTUs,onlyBeta=T){
  
  Betas=lapply(similarity_slices,getBDTT,tree=tree,sampleOTUs=sampleOTUs,onlyBeta=onlyBeta)
  names(Betas)=similarity_slices
  res=do.call(function(...){abind(...,along=0)},Betas)
  return(res)
}

## Examples

library(ape)
library(castor)
library(abind)
library(Matrix)

hist(get_all_node_depths(phy_tree(ps2)))

slices=seq(0,1.6,0.1)

Betas=BDTT(similarity_slices = slices, tree = phy_tree(ps2), sampleOTUs = t(otu_table(ps2)))


# Construct sequence table

predictors="psyllid_spp"
StatsRes=expand.grid(similarity_slices=as.character(slices),predictors=predictors,metric=c("Jac","Bray","Aitch"))
StatsRes[["F.Model"]]=StatsRes[["R2"]]=StatsRes[["Pr(>F)"]]=NA
head(StatsRes)

## ADONIS
for (i in as.character(slices)){
  res=unlist(adonis(Betas[i,"Jac",samples,samples]~ psyllid_spp, method = "euclidean",
                           data = metadata)$aov.tab[1,c(6,5,4)])
  StatsRes[(StatsRes$metric=="Jac")&(StatsRes$similarity_slices==i),4:6]=res
  res=unlist(adonis(Betas[i,"Bray",samples,samples]~ psyllid_spp, method = "euclidean",
                    data = metadata)$aov.tab[1,c(6,5,4)])
  StatsRes[(StatsRes$metric=="Bray")&(StatsRes$similarity_slices==i),4:6]=res
  res=unlist(adonis(Betas[i,"Aitch",samples,samples]~ psyllid_spp, method = "euclidean",
                    data = metadata)$aov.tab[1,c(6,5,4)])
  StatsRes[(StatsRes$metric=="Aitch")&(StatsRes$similarity_slices==i),4:6]=res
}

# Profile for all
ggplot(aes(y=R2,x=similarity_slices,colour=predictors,group=factor(predictors)),data=StatsRes)+geom_point()+geom_line()+facet_wrap(~metric)

## MANTEL TESTS

dists <- expand_grid(dist1 = c("phylo.dist","plant.dist", "spat.dist"), dist2=c("phylo.dist","plant.dist", "spat.dist"))
mantlist_distances <- vector("list", length=length(dists))

for (d in 1:length(dists)){
  test.dist <- get(dists[d])
  
  samples_to_use <- samples[samples %in% colnames(as.matrix(test.dist))]
  
  manlist <- vector("list", length=length(slices))
  for (i in 1:length(slices)){
    
    nslice <- as.character(slices[i])
    #Jac
    mantel <- mantel(Betas[nslice,"Jac",samples_to_use,samples_to_use], test.dist)
    man_jac <- data.frame(metric="Jac", test=dists[d], stat=mantel$statistic,  signif=mantel$signif)
    #Bray
    mantel <- mantel(Betas[nslice,"Bray",samples_to_use,samples_to_use], test.dist)
    man_bray <- data.frame(metric="Bray", test=dists[d], stat=mantel$statistic,  signif=mantel$signif)
    #Aitch_psedo
    mantel <- mantel(Betas[nslice,"Aitch_pseudo",samples_to_use,samples_to_use], test.dist)
    man_aitch_pseudo <- data.frame(metric="Aitch_pseudo", test=dists[d], stat=mantel$statistic,  signif=mantel$signif)
    #Aitch_rep
    mantel <- mantel(Betas[nslice,"Aitch_rep",samples_to_use,samples_to_use], test.dist)
    man_aitch_rep <- data.frame(metric="Aitch_rep", test=dists[d], stat=mantel$statistic,  signif=mantel$signif)
    manlist[[i]] <- bind_rows(man_jac, man_bray, man_aitch_pseudo, man_aitch_rep) %>%
      mutate(slice=slices[i])
  }
  mantlist_distances[[d]] <- bind_rows(manlist)
}

matrix_stats <- bind_rows(mantlist_distances) %>%
  mutate(signif_test = case_when(
    signif < 0.01 ~ TRUE,
    signif > 0.01 ~ FALSE
  ))
gg.mantel <- ggplot(matrix_stats, aes(y=stat,x=slice, colour=test, group=factor(test)))+
  geom_point(aes(shape=signif_test)) +
  geom_line()+ 
  geom_hline(yintercept=0) +
  facet_wrap(~metric) 

## Partial mantels
dists <- expand_grid(dist1 = c("phylo.dist","plant.dist", "spat.dist"), dist2=c("phylo.dist","plant.dist", "spat.dist")) %>%
  filter(!dist1==dist2)
mantlist_distances <- vector("list", length=length(dists))

for (d in 1:nrow(dists)){
  test.dist <- get(dists$dist1[d])
  control.dist <- get(dists$dist2[d])
  samples_to_use <- samples[samples %in% colnames(as.matrix(test.dist))]
  
  manlist <- vector("list", length=length(slices))
  for (i in 1:length(slices)){
    
    nslice <- as.character(slices[i])
    #Jac
    PM <- mantel.partial(Betas[nslice,"Jac",samples_to_use,samples_to_use], test.dist, control.dist)
    man_jac <- data.frame(metric="Jac", dist1=dists$dist1[d], dist2=dists$dist2[d], stat=PM$statistic,  signif=PM$signif)
    #Bray
    PM <- mantel.partial(Betas[nslice,"Bray",samples_to_use,samples_to_use], test.dist, control.dist)
    man_bray <- data.frame(metric="Bray", dist1=dists$dist1[d], dist2=dists$dist2[d], stat=PM$statistic,  signif=PM$signif)
    #Aitch_psedo
    PM <- mantel.partial(Betas[nslice,"Aitch_pseudo",samples_to_use,samples_to_use], test.dist, control.dist)
    man_aitch_pseudo <- data.frame(metric="Aitch_pseudo", dist1=dists$dist1[d], dist2=dists$dist2[d], stat=PM$statistic,  signif=PM$signif)
    #Aitch_rep
    PM <- mantel.partial(Betas[nslice,"Aitch_rep",samples_to_use,samples_to_use], test.dist, control.dist)
    man_aitch_rep <- data.frame(metric="Aitch_rep", dist1=dists$dist1[d], dist2=dists$dist2[d], stat=PM$statistic,  signif=PM$signif)
    manlist[[i]] <- bind_rows(man_jac, man_bray, man_aitch_pseudo, man_aitch_rep) %>%
      mutate(slice=slices[i])
  }
  mantlist_distances[[d]] <- bind_rows(manlist)
}

matrix_stats <- bind_rows(mantlist_distances) %>%
  mutate(signif_test = case_when(
    signif < 0.01 ~ TRUE,
    signif > 0.01 ~ FALSE
  )) %>%
  mutate(test = paste0("microbe~",dist1,"~",dist2))
gg.mantel <- ggplot(matrix_stats, aes(y=stat,x=slice, colour=test, group=factor(test)))+
  geom_point(aes(shape=signif_test)) +
  geom_line()+ 
  geom_hline(yintercept=0) +
  facet_wrap(~metric) 
## Can make all possible combinations using expand grid

# Compare trees

treelist_castor <- vector("list", length=length(slices))
for (i in 1:length(slices)){
  newtree <- collapse_tree_at_resolution(phy_tree(ps2),slices[i], shorten = TRUE)$tree
  treelist_castor[[i]] <- newtree
}

class(treelist_castor) <- "multiPhylo"
gg.castortree <- ggtree(treelist_castor) + facet_wrap(~.id, scale="free") + theme_void()

gg.castortree / gg.mantel

treelist_phyglom <- vector("list", length=length(slices))
for (i in 1:length(slices)){
  newtree <- phy_glom(ps2, similarity=slices[i])
  treelist_phyglom[[i]] <- newtree
}
  
treelist_phyglom <- lapply(treelist_phyglom, phy_tree)

class(treelist_phyglom) <- "multiPhylo"
gg.phyglomtree <- ggtree(treelist_phyglom) + facet_wrap(~.id, scale="free", nrow=1) + theme_void()

gg.castortree / gg.phyglomtree


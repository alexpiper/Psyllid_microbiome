# Copyright (C) 2020 Alexander M Piper
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: alexander.piper@agriculture.vic.gov.au

# Beta diversity through time ---------------------------------------------

BDTT <- function(ps, slices, metrics=c("Bray","Jaccard","Aitchison","Philr"), zeroes="Impute",  multithread = FALSE, quiet=FALSE, progress=FALSE){
  #Check inputs:
  if (is.null(phy_tree(ps))) {stop("A phyloseq object with a phylogenetic tree in the phy_tree slot is required")}
  if(!all(metrics %in% c("Bray", "Jaccard", "Aitchison", "Philr", "Unifrac", "WUnifrac"))){
    stop("metrics must be one or more of: 'Bray', 'Jaccard', 'Aitchison', 'Philr', 'Unifrac' or 'Wunifrac'")}

  #Define getBDTT function
  getBDTT <- function(resolution, ps, metrics=c("Bray","Jaccard","Aitchison","Philr"), zeroes="Impute", parallel=FALSE, quiet=FALSE){
    ps_new <- speedyseq::tree_glom(ps, resolution)
    phy_tree(ps_new) <- multi2di(phy_tree(ps_new))
    phy_tree(ps_new) <- makeNodeLabel(phy_tree(ps_new), method="number", prefix='n')
    message(ntaxa(ps_new), " otus created from a phylogenetic slice at ", resolution)
    if(taxa_are_rows(ps_new)){
      otutab=t(otu_table(ps_new))
    } else if(!taxa_are_rows(ps_new)){
      otutab=otu_table(ps_new)
    }
    # If compositional stats - Deal with zeroes
    if (any(c("Aitchison","Philr") %in% metrics) & stringr::str_to_lower(zeroes) =="pseudo"){
      if(!quiet) {message("Adding a pseudocount for compositional distances")}
      otutab_n0 <- prop.table(otutab + 1, margin=1)
    } else if (any(c("Aitchison","Philr") %in% metrics) & stringr::str_to_lower(zeroes) =="impute"){
      if(!quiet) {message("Imputing zeroes for compositional distances")}
      otutab_n0 <- as.matrix(zCompositions::cmultRepl(otutab, method="CZM", output="prop", suppress.print=quiet))
    }
    ## Calculate distances
    # Bray curtis distance
    if ("Bray" %in% metrics){
      Bray <- as.matrix(vegdist(otutab, method="bray"))
    }
    #Jaccard distance
    if ("Jaccard" %in% metrics){
      Jaccard <- as.matrix(vegdist(otutab, method="jac",binary = T))
    }
    # Aitchison distance
    if ("Aitchison" %in% metrics){
      Aitchison <- as.matrix(vegdist(CoDaSeq::codaSeq.clr(otutab_n0), method="euclidean")) # Replace this with custom CLR function
    }
    #PhilR distance 
    if ("Philr" %in% metrics){
      if(!quiet){
        message("Calculating PhilR distance")
        Philr <- as.matrix(vegdist(philr::philr(otutab_n0, phy_tree(ps_new),
                                                part.weights='enorm.x.gm.counts',
                                                ilr.weights='blw.sqrt'), method="euclidean", na.rm=TRUE))
        message("PhilR distance complete")
      } else if(quiet){
        suppressMessages(
          Philr <- as.matrix(vegdist(philr::philr(otutab_n0, phy_tree(ps_new),
                                                           part.weights='enorm.x.gm.counts',
                                                           ilr.weights='blw.sqrt'), method="euclidean", na.rm=TRUE))
        )
      }
    }
    # Unweighted Unifrac distance
    if ("Unifrac" %in% metrics){
      if(!quiet) {message("Calculating Unweighted Unifrac distance...")}
      Unifrac <- as.matrix(phyloseq::UniFrac(ps_new, weighted=FALSE, parallel = parallel))
      if(!quiet) {message("Unweighted Unifrac distance complete")}
    }
    # Weighted Unifrac distance
    if ("WUnifrac" %in% metrics){
      if(!quiet) {message("Calculating Weighted Unifrac distance...")}
      WUnifrac <- as.matrix(phyloseq::UniFrac(ps_new, weighted=TRUE, parallel = parallel))
      if(!quiet) {message("Unweighted Unifrac distance complete")}
    }
    # Abind only those that are in metrics - Do i have to abind?
    AllBetas <- lapply(metrics, get, envir=sys.frame(sys.parent(0)))
    names(AllBetas) <- metrics
    out <- AllBetas
    return(out)
  }
  
  # setup multithreading
  ncores <- future::availableCores() -1
  if(isTRUE(multithread)){
    cores <- ncores
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multiprocess, workers=cores)
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      message("Warning: the value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multiprocess, workers=cores)
  } else if(isFALSE(multithread) | multithread==1){
    future::plan(future::sequential)
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )
  
  #Apply filt_phmm to all sequences
  Betas <-  furrr::future_map(slices, getBDTT, ps=ps, metrics=metrics, zeroes=zeroes, parallel=FALSE, quiet=quiet, .progress = progress)
  
  #Close all worker threads
  future::plan(future::sequential)

  names(Betas) <- slices
  out <- Betas
  return(out)
}


# Sloans neutral model through time ---------------------------------------

NMTT <- function(ps, slices, return="summary", quiet=FALSE){
  if (is.null(phy_tree(ps))) {stop("A phyloseq object with a phylogenetic tree in the phy_tree slot is required")}
  
  getNMTT <- function(resolution, ps, quiet=FALSE){
    ps_new <- speedyseq::tree_glom(ps, resolution)
    message(ntaxa(ps_new), " otus created from a phylogenetic slice at ", resolution)
    if(taxa_are_rows(ps_new)){
      otutab=t(otu_table(ps_new))
    } else if(!taxa_are_rows(ps_new)){
      otutab=otu_table(ps_new)
    }
    out <- reltools::fit_sncm(otutab, pool=otutab)
    return(out)
  }
  
  neutral <- lapply(slices, getNMTT, ps=ps, quiet=quiet)
  names(neutral) <-slices
  
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

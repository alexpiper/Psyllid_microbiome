# Get Plant Tree ----------------------------------------------------------

library(tidyverse)
library(brranching)
samdf <- read_csv("sample_data/Sample_info2.csv") %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  magrittr::set_rownames(.$SampleID)

plants <- samdf  %>%
  dplyr::select(hostplant_spp)%>%
  unique() %>%
  mutate(branching = brranching::phylomatic_names(.$hostplant_spp, 
                                                format = "isubmit", 
                                                db = "ncbi"))

# Get tree
plant.tree <- brranching::phylomatic(plants$branching, taxnames = TRUE)

# add branch lengths using the wikstrom ages file
plant.tree <- brranching::rbladj(plant.tree, "wikstrom.ages")
plant.tree$tip.label <- str_to_sentence(plant.tree$tip.label)
write.tree(plant.tree, file="sample_data/plant_tree.nwk")
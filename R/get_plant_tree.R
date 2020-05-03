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
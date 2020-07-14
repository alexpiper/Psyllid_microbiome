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

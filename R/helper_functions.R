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
    pull(values) %>%
    mean(na.rm=TRUE)
}

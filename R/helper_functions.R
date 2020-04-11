# Negative to positive eig ------------------------------------------------
neg_to_pos_eig = function(PCo){
  min_val = PCo %>% as.vector %>% min
  print(min_val)
  if(min_val < 0){
    PCo = PCo + abs(min_val)
  }
  return(PCo)
}

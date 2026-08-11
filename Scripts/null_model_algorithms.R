#------------------------------------------------------------------------------
#
#  Null model algorithm functions
#
#------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 2-13-2025

# Description: A set of functions to randomize community, trait, or phylogenetic
# data for null model creation

# Alpha diversity null models --------------------------------------------------

# Region taxa swap
taxa_swap_region <- function(com,region_bridge) {
  
  # All species list
  sp <- colnames(com)[order(colnames(com))]
  comid <- row.names(com)
  
  # Region list
  regions <- region_bridge %>% distinct(HUC_2) %>% pull()
  
  # Create null models for each region
  shuffle_list <- lapply(regions, function(x){
    
    # Create regional species pool
    reg_com <- com %>% 
      tibble::rownames_to_column("COMID") %>% 
      left_join(region_bridge,by = join_by(COMID)) %>% 
      filter(HUC_2 == x) %>% 
      select(-HUC_12,-HUC_2,-HUC_8) %>% 
      select(where(~ any(. != 0))) %>% 
      tibble::column_to_rownames("COMID")
    sp_pool <- colnames(reg_com)
    
    # Shuffle regional taxa labels
    null_sp <- sample(sp_pool,length(sp_pool))
    colnames(reg_com) <- null_sp
    
    # Create absense (0) data for species outside of pool
    abs <- setdiff(sp,sp_pool)
    reg_com[,abs] <- 0
    stopifnot(ncol(reg_com) ==ncol(com))
    
    # Reorder for rowbinding
    reg_com[,sp]
  })
  
  # bind regional null models into one dataset
  bind_rows(shuffle_list)[comid,]
  
  
}

# Beta diversity null models  --------------------------------------------------

# Functional null model
trait_swap <- function(trait) {
  sp <- row.names(trait)
  n_sp <- length(sp)
  sp_null <-sample(sp,n_sp)
  row.names(trait) <- sp_null
  trait[order(sp_null),]
}

# Phylogenetic null model
tree_swap <- function(tree) {
  sp <- tree$tip.label
  n_sp <- length(sp)
  sp_null <-sample(sp,n_sp)
  tree$tip.label <- sp_null
  tree
}



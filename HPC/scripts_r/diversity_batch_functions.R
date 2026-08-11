#-------------------------------------------------------------------------------
#
#   Diversity batch functions
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 03/16/2026

# Description: Functions to batch estimate functional and phylogenetic diversity
# metrics for the use in null models.


# Beta functions  --------------------------------------------------------------
#' Batch estimate kernel density based  functional beta diversity and lcbd
#' 
#' @description Estimates functional beta diversity and local contribution to
#' beta diversity (LCBD) using kernel density hypervolumes. Designed to run
#' within the parallel function
#' 
#' @param com Community matrix with columns for species and rows for site
#' @param trait_hyp Trait hypervolume created using BAT::hyper.build 
#' @inheritParams BAT::kernel.build
#' 
#' @returns Named list containing pairwise beta diversity output list from
#' BAT::kernel.beta and a data.frame containing LCBD values for total, 
#' replacement, and richness difference components
#' 

kernel_beta_batch <- function(com,trait_hyp,cores =1, abund = F){
  
  ### Check that species are in correct order
  match <-trait_match(com,trait_hyp)
  com_match <- match[[1]]
  trait_match <- match[[2]]
  
  ### Create hyper volumes ###
  kernel <- BAT::kernel.build(com_match,
                              trait_match,
                              abund = abund,
                              cores = cores)
  
  ### Calculate beta diversity  ###
  beta <- BAT::kernel.beta(kernel,func = "sorensen")[1:3]
  
  # Create data frame with LCBD of all three components
  lcbd <- purrr::map2(beta,names(beta),lcbd_batch, dim = "fun") %>% 
    purrr::reduce(left_join,by = join_by(COMID)) %>% 
    tibble::column_to_rownames("COMID")
  
  list(beta,lcbd)
  
}


#' Batch estimate taxonomic beta diversity and lcbd
#' 
#' @description Estimates taxonomic beta diversity and local contribution to
#' beta diversity (LCBD). Designed to run within the parallel function
#' 
#' @param com Community matrix with columns for species and rows for site
#' @inheritParams BAT::beta
#' 
#' @returns Named list containing pairwise beta diversity output list from
#' BAT::beta and a data.frame containing LCBD values for total, 
#' replacement, and richness difference components
#' 

tax_beta_batch <- function(com,abund = F){
  
  ### Calculate beta diversity  ###
  beta <- BAT::beta(as.matrix(com),func = "sorensen",abund = abund)[1:3]
  
  # Create data frame with LCBD of all three components
  lcbd <- purrr::map2(beta,names(beta),lcbd_batch,dim="tax") %>% 
    purrr::reduce(left_join,by = join_by(COMID)) %>% 
    tibble::column_to_rownames("COMID")
  
  list(beta,lcbd)
  
}


#' Batch estimate phylogenetic beta diversity and lcbd
#' 
#' @description Estimates phylogenetic beta diversity and local contribution to
#' beta diversity (LCBD). Designed to run within the parallel function
#' 
#' @param com Community matrix with columns for species and rows for site
#' @inheritParams BAT::beta
#' 
#' @returns Named list containing pairwise beta diversity output list from
#' BAT::beta and a data.frame containing LCBD values for total, 
#' replacement, and richness difference components
#' 

phy_beta_batch <- function(com,tree,abund = F){
  
  ### Calculate beta diversity  ###
  beta <- BAT::beta(com,tree,func = "sorensen",abund = abund)[1:3]
  
  # Create data frame with LCBD of all three components
  lcbd <- purrr::map2(beta,names(beta),lcbd_batch,dim="phy") %>% 
    purrr::reduce(left_join,by = join_by(COMID)) %>% 
    tibble::column_to_rownames("COMID")
  
  list(beta,lcbd)
  
}


# Alpha functions  -------------------------------------------------------------

#' Batch estimate kernel density based  functional alpha diversity
#' 
#' @description Estimates functional alpha diversity using kernel density 
#' hypervolumes. Designed to run within the parallel function.
#' 
#' @param com Community matrix with columns for species and rows for site
#' @param trait_hyp Trait hypervolume created using BAT::hyper.build 
#' @inheritParams BAT::kernel.build
#' 
#' @returns named vector containing diversity values
#' 

kernel_alpha_batch <- function(com,trait_hyp,cores =1, abund = F){
  
  ### Check that species are in correct order
  match <-trait_match(com,trait_hyp)
  com_match <- match[[1]]
  trait_match <- match[[2]]
  
  ### Create hyper volumes ###
  kernel <- BAT::kernel.build(com_match,
                              trait_match,
                              abund = abund,
                              cores = cores)
  
  ### Calculate alpha diversity  ###
  BAT::kernel.alpha(kernel)
  
}


# helper functions  ------------------------------------------------------------

# Batch convert LCBD output to dataframe
lcbd_batch <- function(D,dim,name) {
  
  # check if euclidean
  sqrt <- ifelse(ade4::is.euclid(D) == T, F, T)
  
  # estimate LCBD
  lcbd <- adespatial::LCBD.comp(D,sqrt.D = sqrt)$LCBD
  
  # Create output dataframe
  lcbd_colname <- paste0(dim,"_",name)
  out <- data.frame(COMID = labels(D))
  out[,lcbd_colname] <- lcbd
  row.names(out) <- out$COMID
  out
}

# Reorder community and trait matrices to match
trait_match <- function(com,trait_hyp){
  sp <- colnames(com)[order(colnames(com))]
  com_match <- com[,sp]
  trait_match <- trait_hyp[sp,]
  
  stopifnot("Species in com do not match species in trait" = 
              all(colnames(com_match) == row.names(trait_match)))
  
  list(com_match,trait_match)
}


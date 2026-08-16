#' Regional constrained taxa-swap null model algorithm
#' 
#' @description Random shuffles species labels in community matrix, constraining
#' swaps to some defined regional species pool (HUC_2)
#' 
#' @param com Community matrix with columns for species and rows for site
#' @param region_bridge  data.frame with column COMID that matches row names
#' of com, and HUC_2 which define that site's region.
#'
#' @returns randomly shuffled community matrix with identical labels and 
#' dimensions to input matrix
#' 

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


# Unconstrained taxa swap helpers  ---------------------------------------------
.comm_swap <-function(comm) {
  sp <- row.names(trait)
  n_sp <- length(sp)
  sp_null <-sample(sp,n_sp)
  row.names(trait) <- sp_null
  trait[order(sp_null),]
}

.trait_swap <- function(trait) {
  sp <- row.names(trait)
  n_sp <- length(sp)
  sp_null <-sample(sp,n_sp)
  row.names(trait) <- sp_null
  trait[order(sp_null),]
}

.tree_swap <- function(tree) {
  sp <- tree$tip.label
  n_sp <- length(sp)
  sp_null <-sample(sp,n_sp)
  tree$tip.label <- sp_null
  tree
}

#' # Tax Beta diversity null models  ----------------------------------------------
#' 
#' #' Nonnative species swap taxonomic beta diversity change null model algorithm
#' #' 
#' #' @description Randomizes the occurrences of nonnative species in a 
#' #' contemporary community data set, for use in taxonomic beta diversity change
#' #' null models.
#' #' 
#' #' @param mod_com Community matrix containing all species (native+nonnative)  
#' #' with columns for species and rows for site. Must have the same row names as
#' #' mod_com.
#' #' @param his_com Community matrix containing only native species with 
#' #' columns for species and rows for site. Must have the same row names as
#' #' mod_com. "Pseudo-historical" data
#' #' @param region_bridge data.frame with column COMID that matches row names
#' #' of com, and HUC_8 which define that site's region, and HUC_8 which defines
#' #' a site's subbasin (used for native vs nonnative identification)
#' #' @param const character ("reg.all","reg.cont", or "reg.freq"). Defines 
#' #' constraints to null model randomization. See details for more information 
#' #' 
#' #'   
#' #' @details 
#' #' Null model algorithm based on Leprieur et al., (2008), that randomizes the 
#' #' occurrences of nonnative species to test if directional  beta diversity 
#' #' change is more extreme then expected by random distributions of nonnative
#' #' species. This algorithm should only be used in pseudo-historic approaches, 
#' #' where full community data is compared to community data with nonnative species 
#' #' withheld as pseudo-historical community data.
#' #' 
#' #' This null model maintains an equiprobable total of columns by randomizing the 
#' #' spatial distribution of nonnative species, maintaining their occurrence 
#' #' frequency, but allowing communities to be unconstrained in the number of
#' #' nonnative species they receive. Leprieur et al., (2008) provides detailed 
#' #' arguments for the ecological validity of this equiprobable null model in 
#' #' the context of invasion-driven Δβ. In short, the authors argue that the 
#' #' model fits the expectation that nonnative species differ in their 
#' #' colonization ability and propagule pressure, most communities are 
#' #' susceptible to invasions, and that communities are rarely saturated. 
#' #' 
#' #' While this algorithm contains the occurrence frequency of each
#' #' nonnative species to observed frequencies, this function allows 
#' #' users to specify additional constraints. "reg.all" allows for a nonnative
#' #' species to be placed into any region (HUC_2) and does not maintain regional
#' #' occurrence frequencies; "reg.cont" restricts nonnative species into regions 
#' #' where they have been observed, nonnative but does not maintain regional
#' #' occurrence frequencies; and "reg.freq" restrict nonnative species
#' #' to regions where they have occurred and maintains regional
#' #' occurrence frequencies.
#' #' 
#' #' @returns community matrix with identical dimensions and labels as mod_com,
#' #' but random distribution of nonnative species occurrences
#' #' 
#' 
#' inv_swap <- function(mod_com,his_com,region_bridge,const = "reg.all") {
#'   
#'   ## Prep community data ##
#'   
#'   # all species in analysis
#'   all_sp <- colnames(mod_com)[order(colnames(mod_com))]
#'   comid <- row.names(mod_com)
#'   
#'   # Format community matrices to have the same species, with zeros for missing
#'   # species
#'   miss_sp <- setdiff(all_sp,colnames(his_com))
#'   his_com_for <- his_com
#'   his_com_for[,miss_sp] <-0
#'   
#'   # Arranage species into same order 
#'   his_com_for <- his_com_for[,all_sp]
#'   mod_com_for <- mod_com[,all_sp]
#'   
#'   # CHeck that species are in same order
#'   stopifnot(
#'     "Species reording failed" = 
#'       all(colnames(mod_com_for) == colnames(his_com_for))
#'   )
#'   
#'   # Create community matrix for only nonnative species
#'   delta_com <- mod_com_for-his_com_for
#'   
#'   # Number of nonnative occurrences for each species
#'   nn_oc <- colSums(delta_com)[colSums(delta_com)>0]
#'   
#'   ## Set null model constraints ##
#'   stopifnot('Const needs to be "reg.all", "reg.cont", or "reg.freq"' =
#'               const %in% c("reg.all","reg.cont","reg.freq"))
#'   
#'   # Not constrained by region, each nn sp can occur in any huc2
#'   if(const == "reg.all") {
#'     
#'     # assign each nonnative species all possible huc2s
#'     all_huc2 <- region_bridge %>% 
#'       distinct(HUC_2) %>% 
#'       pull()
#'     nn_region_list <- lapply(nn_oc, function(x) all_huc2)
#'   }
#'   
#'   # Constrained by region, nn sp only occur in huc2s in which they have been 
#'   # observed. Frequency within each region is not maintained
#'   if(const == "reg.cont") {
#'     
#'     # assign each nn species the huc2 which they have been observed in data
#'     nn_region_df <- delta_com %>% 
#'       select(where(~sum(.) > 0)) %>% 
#'       tibble::rownames_to_column("COMID") %>% 
#'       left_join(region_bridge %>% select(COMID,HUC_2)) %>% 
#'       select(-COMID) %>% 
#'       group_by(HUC_2) %>% 
#'       summarise(across(everything(), sum)) %>% 
#'       tibble::column_to_rownames("HUC_2")
#'     nn_region_list <- lapply(
#'       as.list(nn_region_df), 
#'       function(x) row.names(nn_region_df[x>0,])
#'     )
#'   }
#'   
#'   #Constrained by region, nn sp only occur in huc2s in which they have been 
#'   # observed. Frequency within each region is maintained
#'   if (const == "reg.freq"){
#'     
#'     # Extract the region-specific, nn occurance frequencies 
#'     region_nn_df <- delta_com %>% 
#'       select(where(~sum(.) > 0)) %>% 
#'       tibble::rownames_to_column("COMID") %>% 
#'       left_join(region_bridge %>% select(COMID,HUC_2)) %>% 
#'       select(-COMID) %>% 
#'       group_by(HUC_2) %>% 
#'       summarise(across(everything(), sum)) %>% 
#'       tibble::column_to_rownames("HUC_2") %>% 
#'       t() %>% 
#'       as.data.frame() 
#'     region_nn_oc <-lapply(as.list(region_nn_df),function(x) {
#'       names(x) <- row.names(region_nn_df)
#'       x[x>0,drop =F]
#'     })
#'   }else{
#'     
#'     # Place total occurance frequncies inot list of 1
#'     region_nn_oc <- list(nn_oc)
#'     names(region_nn_oc) <- "all"
#'   }
#'   
#'   
#'   ## Assemble null nonnative communities  ##
#'   nn_com_region <- purrr::map2(
#'     region_nn_oc,names(region_nn_oc), 
#'     function(nn.oc,region) {
#'       
#'       # Extract selected region, if regional occurances are maintained
#'       if(region != "all"){
#'         nn_region_list <- lapply(nn.oc, function(x) region)
#'       } 
#'       
#'       # Create a nonnative species community matrix with random occurrences, 
#'       # maintaining overall occupancy frequencies
#'       nn_com_list <-purrr::map2(nn.oc,names(nn.oc), function(oc,nn_sp){
#'         
#'         # Remove native range from comid list, and limit to desired regional 
#'         # constraints
#'         native_hucs <- native_range[[nn_sp]]
#'         native_comids <- row.names(
#'           his_com_for[his_com_for[,nn_sp] ==1,nn_sp,drop = F]
#'         ) 
#'         nn_hucs <- nn_region_list[[nn_sp]]
#'         nn_comids <- region_bridge %>% 
#'           filter(!HUC_8 %in% native_hucs,
#'                  !COMID %in% native_comids,
#'                  HUC_2 %in% nn_hucs) %>% 
#'           pull(COMID)
#'         
#'         # Randomly assign species to nonnative comids based on occ. freq
#'         nn_comid <- sample(nn_comids,oc)
#'         nn_df <- data.frame(COMID = nn_comid)
#'         nn_df[,nn_sp] <- 1 
#'         nn_df
#'       })
#'       
#'       # Merge all nonnative species into data frame
#'       nn_com <- purrr::reduce(nn_com_list,full_join,by = join_by(COMID)) %>% 
#'         column_to_rownames("COMID")  
#'       nn_com[is.na(nn_com)] <- 0  # change NAs to absences (0)
#'       nn_com
#'     })
#'   
#'   # Merge regional data frames, if regional occrences maintained
#'   nn_com <- bind_rows(nn_com_region)
#'   nn_com[is.na(nn_com)] <- 0
#'   
#'   ## Reassemble null contemporary communities  ##
#'   
#'   # Create a list of each nonnative species at a comid
#'   nn_com_t <- as.data.frame(t(nn_com))
#'   nn_sp_list <- lapply(
#'     as.list(nn_com_t), 
#'     function(x) row.names(nn_com_t[x==1,])
#'   )
#'   
#'   # For loop to add nonnative species to each invaded comid (for loop in this 
#'   #case is more parsimonous then lapply)
#'   null_com <- his_com_for
#'   for(i in 1:length(nn_sp_list)) {
#'     x <- names(nn_sp_list)[i]
#'     y <- nn_sp_list[[i]]
#'     stopifnot("nn species assigned to native range" = all(null_com[x,y] ==0))
#'     null_com[x,y] <- 1
#'   }
#'   null_com
#' }
#' 




#-------------------------------------------------------------------------------
#
#  Null model input data creation
#
#-------------------------------------------------------------------------------

# Author: 

# Created: 2-13-2025

# Description: Prepares input data for the calculation of observed and null
# diversity values on high performance computer clusters. 


# House Keeping ----------------------------------------------------------------
rm(list = ls())

# Packages
library(BAT)
library(dplyr)
library(tibble)
library(purrr)
library(ape)

# data directories
data_dir <- "Diversity Input Data"

# Export directories
beta_dir <- "HPC/beta_null_input_data"
obs_dir <- "HPC/beta_obs_input_data"
alpha_dir <- "HPC/alpha_null_input_data"


# load in custom null model algorithms
algo_dir <-"HPC/scripts_do_not_run"
source(file.path(algo_dir,"null_model_algorithms.R"))

# Load data
mod_com <- readRDS(file.path(data_dir,"mod_com_diversity_input.rds"))
his_com <- readRDS(file.path(data_dir,"his_com_diversity_input.rds"))
traits <- readRDS(file.path(data_dir,"trait_diversity_input.rds"))
tree <- readRDS(file.path(data_dir,"phylo_tree.rds"))
max_iter  <- 999


# Taxonomic and phylo data prep ------------------------------------------------

# Check for sites with no species. If one site has no species this will still 
# allow total beta diversity to be calculated but prevent the calculation of 
# components. If more than one site has no species than this will cause NAs 
# for Beta diversity between those sites!
nrow(mod_com[rowSums(mod_com[,-1]) == 0,])
nrow(his_com[rowSums(his_com[,-1]) == 0,])

# Create a COMID to HUC12 bridge df. For regional species pools. Skip this code
# if regional species pools are not being used. Modify to whatever rownames
# and region names used
region_bridge <- mod_com %>% 
  select(HUC_12) %>% 
  tibble::rownames_to_column("COMID") %>% 
  mutate(HUC_2 = substr(HUC_12,1,2),
         HUC_8 = substr(HUC_12,1,8))
mod_com <- mod_com %>% select(-HUC_12)
his_com <- his_com %>% select(-HUC_12)

#  Make sure row names are what expected
row.names(mod_com)
row.names(his_com)


### Phylogenetic ###
# Remove all species not in analysis
sp_list <- unique(colnames(mod_com),colnames(his_com))
tree <- ape::keep.tip(tree,sp_list)


# Functional data prep ---------------------------------------------------------

# removal all species not included in full community data
traits <- traits[row.names(traits) %in% sp_list,] 

# Final check to ensure all species are in trait data
sp_list[!sp_list %in% row.names(traits)]

# Reduce dimensionality and create hypervolume (use code below to determine 
# number of axes to retain)
trait_gower <- BAT::gower(traits)
trait_hyp <- BAT::hyper.build(trait_gower,axes = 4)  # reduce using pcoa 

# If using PCoA to reduce dimensionality, determine number of axes to retain
# using broken stick
#
# # calculate gower distance of all species
# trait_dist <- BAT::gower(traits,"gower")
# 
# 
# # conduct pcoa
# trait_pca <- ape::pcoa(trait_dist)
# # Export eigen table
# write.csv(trait_pca[["values"]],"syndromee_eigen_table.csv")
# 
# ## Select number of pc axes with broken stick ##
# # Extract eigenvalues from the PCA result
# eigenvalues <- trait_pca$values$Relative_eig
# 
# # Function to calculate broken stick values
# broken_stick <- function(eigenvalues) {
#   n <- length(eigenvalues)
#   bs_values <- sapply(1:n, function(i) sum(1/i:n))
#   bs_values <- bs_values / sum(bs_values)
#   return(bs_values)
# }
# 
# # Calculate broken stick values
# bs_values <- broken_stick(eigenvalues)
#
# # Plot eigenvalues and broken stick values for comparison
# 
# tiff(filename = "syndrome_trait_bs.tiff",width = 1600,height = 1200)
# par(mar=c(8,8,4,2)+.8)
# par(mgp=c(6,2,0))
# plot(eigenvalues,
#      type = "b",
#      col = "blue",
#      xlab = "Principal Components",
#      ylab = "% Eigenvalues",
#      cex.lab = 3,
#      cex.axis = 2,
#      cex = 2,
#      lwd = 2,
#      xaxt='n')
# axis(side=1, at=c(0:19),cex.axis = 2)
# lines(bs_values,
#       type = "b",
#       col = "red",
#       cex = 2,
#       lwd = 2)
# 
# # Add legend
# legend("topright",
#        legend = c("% Eigenvalues", "Broken Stick"),
#        col = c("blue", "red"),
#        lty = 1,
#        cex = 3,
#        lwd =2)
# dev.off()
# 
# # Determine the number of axes to retain
# num_axes_to_retain <- sum(eigenvalues > bs_values)
# cat("Number of axes to retain:", num_axes_to_retain, "\n")
# ## 4 axes should be retained according to broken stick
# # 0.6456354 varation explained
#


# Create functional space plots-------------------------------------------------
# library(corrplot)
# 
# ## Data prep ##
# trait_pcs <- as.data.frame(trait_pca$vectors[,1:4]) %>%  # extract 4 pcs
#   tibble::rownames_to_column(var = "Sci_name") %>%  
#   left_join(traits %>%
#               tibble::rownames_to_column(var = "Sci_name")) %>% 
#   
#   # reodrer traits to matach table in appendix
#   select(
#     Axis.1,Axis.2,Axis.3,Axis.4,MAXBODYL,FECUND,LONG,AGEMAT,
#     PARENTC,EggSize,`Herbivore-detritivore`,Omnivore,Invertivore,
#     `Invertivore-piscivore`,Piscivore,Benthic,`Non-benthic`,SLOWCURR,MODCURR,
#          FASTCURR,MAXTEMP,Rubble,Sand,`Silt/mud`,Vegetation,Various
#     ) %>%  
#   rename("PC1" = Axis.1,
#          "PC2" = Axis.2,
#          "PC3" = Axis.3,
#          "PC4" = Axis.4,
#          "Max total length" = MAXBODYL,
#          "Fecundity" = FECUND,
#          "Longevity" = LONG,
#          "Age at maturity" = AGEMAT,
#          "Parental care" = PARENTC,
#          "Egg size" = EggSize,
#          "Slow Current" = SLOWCURR,
#          "Moderate Current" = MODCURR,
#          "Fast Current" = FASTCURR,
#          "Max temperature" = MAXTEMP)  # rename traits for better presentation
# 
# ## Biplots ##
# pc_list <- c("PC1","PC2","PC3","PC4")  # name of each pc axis
# pc_comb <- combn(pc_list,2)  # every combo of pc axes in a dataframe
# 
# For loop to create plots
# for (i in 1:ncol(pc_comb)){
# 
#   png(
#     filename = paste(
#     "trait_pcoa_plots/trait_pcoa_",
#      pc_comb[1,i],"_",
#      pc_comb[2,i],
#      ".png",
#      sep = ""
#      ), 
#     width = 1600, 
#     height = 1600
#     )  
#   par(mar=c(8,8,4,2)+1.2)
#   par(mgp=c(6,5,0))
#   # plot every combination of pcs
#   plot(trait_pcs[,pc_comb[2,i]] ~ trait_pcs[,pc_comb[1,i]],  
#        xlab = "",  # remove x label
#        ylab = "",  # remove ylabel
#        pch = 16,  # point type
#        cex = 5,  # point size
#        cex.axis = 7)  # text size)
#   abline(v = 0, lty = 2, lwd = 5)
#   abline(h = 0, lty =2, lwd = 5)
#   dev.off()
# }
# 
# ## Corplots
# # calcualte pearson correlation between pc axes and traits
# # trait_cors <- cor(trait_pcs)[5:nrow(cor(trait_pcs)),1:4]  
# #
# # # Plot correlations
# # png(filename = "syndrome_trait_cor.tiff",width = 800,height = 1600)
# # corrplot(trait_cors,
# #          addgrid.col= "black",  # change grid to black
# #          tl.col = "black",  # change text to black
# #          tl.cex = 3,  # text size
# #          cl.pos = "r",  # legend size
# #          cl.ratio = .5,  # legend width
# #          cl.cex = 3,  # legend text size
# #          cl.length = 3) # lengend scaling
# # dev.off()
# 


# Observed input data  ---------------------------------------------------------

# Phylogenetic
mod_phy_obs <- list(time="mod",com=mod_com,tree=tree)
his_phy_obs <- list(time="his",com=his_com,tree=tree)

# Functional
mod_fun_obs <- list(time="mod",com=mod_com,trait=trait_hyp)
his_fun_obs <- list(time="his",com=his_com,trait=trait_hyp)


# Null alpha input data --------------------------------------------------------

# Create null native community matrices
alpha_com_list <- parallel::mclapply(
  seq_len(max_iter),
  function(x) taxa_swap_region(his_com,region_bridge),
  mc.cores = 10
  )
names(alpha_com_list) <- sprintf("%03d",seq(1:length(alpha_com_list)))

# Create input lists
phy_alpha_input <- list(
  time="his",
  com_list = alpha_com_list, 
  tree = tree
  )

# Create splits in trait list based on number of nodes needed
n_groups_f <- 500
grouping_f <-cut(
  seq(1:max_iter),
  breaks = n_groups_f,
  labels = F
  )
com_null_split <- split(alpha_com_list,grouping_f)

fun_alpha_input <- lapply(
  com_null_split, 
  function(x) list(time="his",com_list=x,trait=trait_hyp)
  )


# Null beta input data - taxonomic  --------------------------------------------
tax_beta_input <- parallel::mclapply(
  1:max_iter, 
  function(x) inv_swap(
    mod_com,
    his_com,
    region_bridge,
    "reg.cont"
    ),
  mc.cores =12
  )


# Null phylogenetic beta input data --------------------------------------------

# Create null functional hyper volumes
phy_null_list <- replicate(max_iter,tree_swap(tree),simplify = F)
names(phy_null_list) <- sprintf("%03d",seq(1:length(phy_null_list)))

# How many cores per nodes to get 111 total
sapply(c(1, 3,  37, 111),function(x) 111/x)
# Could do 3 cores of 37 cores this will take 9 times the run of one iteration 
# (i.e. ~ 9 hours)

# Create splits in trait list based on number of nodes needed
n_groups_p <- 3
grouping_p <-cut(seq(1:max_iter),breaks = n_groups_p,labels = F)
phy_null_split <- split(phy_null_list,grouping_p)

# Add community data to input lists for both contemporary and native pools
mod_phy_input <- lapply(
  phy_null_split, 
  function(x) list(time="mod",com=mod_com,phy_list=x)
  )
his_phy_input <- lapply(
  phy_null_split, 
  function(x) list(time="his",com=his_com,phy_list=x)
  )

# Merge time periods into one list
names(his_phy_input) <- as.numeric(names(his_phy_input)) +n_groups_p
phy_input <- append(mod_phy_input,his_phy_input)


# Null Functional beta input data ----------------------------------------------

# Create null functional hypervolumes
trait_null_list <- replicate(max_iter,trait_swap(trait_hyp),simplify = F)
names(trait_null_list) <- sprintf("%03d",seq(1:length(trait_null_list)))

# How many cores per nodes to get 1000 total
sapply(c(1, 3, 9, 27, 37, 111),function(x) 999/x)
# Could do 37 nodes with 27cores  or 27 nodes with 37 cores
# Double this for the total number of cores need for both time periods
# will need to re request 52 or 74 nodes total

# Create splits in trait list based on number of nodes needed
n_groups_f <- 500
grouping_f <-cut(seq(1:max_iter),breaks = n_groups_f,labels = F)
trait_null_split <- split(trait_null_list,grouping_f)

# Add community data to input lists for both contemporary and native pools
mod_fun_input <- lapply(trait_null_split, 
                        function(x) list(time="mod",com=mod_com,trait_list=x))
his_fun_input <- lapply(trait_null_split, 
                        function(x) list(time="his",com=his_com,trait_list=x))

# Merge time periods into one list
names(his_fun_input) <- as.numeric(names(his_fun_input)) +n_groups_f
fun_input <- append(mod_fun_input,his_fun_input)


# Export  ----------------------------------------------------------------------

# observed input lists
saveRDS(mod_fun_obs,file.path(obs_dir,"mod_fun_obs_input_list.rds"))
saveRDS(his_fun_obs,file.path(obs_dir,"his_fun_obs_input_list.rds"))
saveRDS(mod_phy_obs,file.path(obs_dir,"mod_phy_obs_input_list.rds"))
saveRDS(his_phy_obs,file.path(obs_dir,"his_phy_obs_input_list.rds"))

# Null alpha input lists
saveRDS(
  phy_alpha_input,
  file.path(alpha_dir,"his_phy_alpha_null_input_list.rds")
  )
purrr::map2(
  fun_alpha_input,names(fun_alpha_input), 
  function(x,y) saveRDS(
    x,
    file.path(alpha_dir,paste0("his_fun_alpha_null_input_list_",y,".rds"))
    )
  )

# Null beta taxa input lists
saveRDS(tax_beta_input,file.path(tax_dir,"tax_null_input_list.rds"))

# Null beta input lists
purrr::map2(
  fun_input,names(fun_input), 
  function(x,y) saveRDS(
    x,
    file.path(beta_dir,paste0("fun_null_input_list_",y,".rds"))
    )
  )
purrr::map2(
  phy_input,
  names(phy_input), 
  function(x,y) saveRDS(
    x,
    file.path(beta_dir,paste0("phy_null_input_list_",y,".rds"))
    )
  )


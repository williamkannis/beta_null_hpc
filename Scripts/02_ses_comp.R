# ------------------------------------------------------------------------------
#
#  SES output compilation  
#
# ------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 04/28/2026

# Description: Compiles and formats the resulting SES, ES, and diagnostic 
# stats across files. Also visualizes normality diagnostics to allow users to 
# decide between SES or ES values for further analyses


# House Keeping ----------------------------------------------------------------
rm(list=ls())

# Directories
out_dir <- "HPC/ses_outputs"
obs_dir <- "HPC/obs_outputs"
export_dir <- "Diversity Output Data"


# Packages
library(dplyr)
library(ggplot2)
library(purrr)


# Load in outputs --------------------------------------------------------------

# all outputs
out_files <- list.files(out_dir)

# Metrics of interst
mets <- c(
  "alpha",  # native alpha
  "his_lcbd",  # native LCBD
  "mod_lcbd",  # contemporary LCBD
  "_d_lcbd",  # change in lcbd
  "his_B",  # native beta diversity
  "mod_B",  # contemporary beta diversity
  "_d_B"  # change in beta diversity
  )

# Extract files for each metric of interest
file_list <- lapply(mets, function (x) out_files[grepl(x,out_files)])
names(file_list) <- mets


# Format outputs  --------------------------------------------------------------
out_list <- lapply(file_list, function(y) {
  
  # Load in files
  ls <- lapply(y,function(x) readRDS(file.path(out_dir,x)))
  names(ls) <- y

  # Transpose list so BLANK is on outside
  ls_t <- purrr::transpose(ls)
  
  # Format lcbd matrices
  if(all(grepl("lcbd",y))) {
    out <- lapply(ls_t, function(x){
      ls <- lapply(x, function(y) {
        as.data.frame(y) %>% tibble::rownames_to_column("COMID")
        }
        )
      purrr::reduce(ls,left_join,by=join_by(COMID)) %>%
        tibble::column_to_rownames("COMID")
    })
    return(out)
    }

  if(all(grepl("alpha",y)|grepl("_B",y))){
    out <- lapply(ls_t, function(x){
      df <-as.data.frame(do.call(cbind,x))
      colnames(df) <- sub("_ses_out.rds","",colnames(df))
      df
    })
    return(out)
  }
  
})

# Observed taxonomic data ------------------------------------------------------

# manually merge in taxonomic observed values

# load in obs tax beta
his_tax_obs <- readRDS(file.path(obs_dir,"his_tax_beta_obs.rds"))
mod_tax_obs <- readRDS(file.path(obs_dir,"mod_tax_beta_obs.rds"))

# Extract beta diversity components into matrix
his_tax_b <- do.call(cbind,lapply(his_tax_obs[[1]],as.vector))
mod_tax_b <- do.call(cbind,lapply(mod_tax_obs[[1]],as.vector))
colnames(his_tax_b) <- paste0("tax_his_",colnames(his_tax_b))
colnames(mod_tax_b) <- paste0("tax_mod_",colnames(mod_tax_b))

# Merge to beta data
out_list$his_B$obs <- cbind(his_tax_b,out_list$his_B$obs)  
out_list$mod_B$obs <- cbind(mod_tax_b,out_list$mod_B$obs)  

# Extract and format LCBD
his_tax_lcbd <- his_tax_obs[[2]]  
mod_tax_lcbd <- mod_tax_obs[[2]]  

# are LCBD rows in same order
all(row.names(his_tax_lcbd) == row.names(mod_tax_lcbd))

# Calculate delta LCBD
d_tax_lcbd <- mod_tax_lcbd - his_tax_lcbd

# Add COMID column to lcbd
his_tax_lcbd <- his_tax_lcbd %>% tibble::rownames_to_column("COMID")
d_tax_lcbd <- d_tax_lcbd %>% tibble::rownames_to_column("COMID")

# Merge to LCBD data
out_list$his_lcbd$obs <- out_list$his_lcbd$obs %>% 
  tibble::rownames_to_column("COMID") %>% 
  left_join(his_tax_lcbd,by = join_by(COMID)) %>% 
  tibble::column_to_rownames("COMID")
out_list$`_d_lcbd`$obs <- out_list$`_d_lcbd`$obs %>% 
  tibble::rownames_to_column("COMID") %>% 
  left_join(d_tax_lcbd,by = join_by(COMID)) %>% 
  tibble::column_to_rownames("COMID")


# Transpose output list
out_list_t <- purrr::transpose(out_list)


# Check for skew in null distributions  ----------------------------------------

# Skew plots
for (i in 1:length(out_list_t$skew)) {
  skew_df <- out_list_t$skew[[i]]
  boxplot(skew_df)
  abline(h=2)
  abline(h=-2)
  abline(h=0)
  title(names(out_list_t$skew)[i])
}


# Summary Stats ----------------------------------------------------------------

# Percent change in observed mean pairwise beta
mpw_pct_list <-lapply(c("mod_B","his_B"), function(i){
  time_name = ifelse(i=="mod_B","Contemporary","Native")
  out_list_t$obs[[i]] %>% 
    summarise(across(everything(),~mean(.x)))
})
do.call(function(x,y) 100*(1-x/y),mpw_pct_list)


# Observed delta LCBD summary
summary(out_list_t$obs$`_d_lcbd`)

# ES delta LCBD summary
summary(out_list_t$empirical_es$`_d_lcbd`)


# Prepare diversity data for analyses  -----------------------------------------

# As an example, we prepare native alpha, native lcbd, and delta lcbd data
# into data.frames for hypothetical analyses.

# Subset to only include variables used in analysis
analysis_list <-out_list[c("alpha","his_lcbd","_d_lcbd")]

final_list <-lapply(analysis_list, function (x) {
  
  # select observed and ES or SES
  sub.vec <- c(
    "obs",
     # "ses",  # if null is not skewed, you can use this instead
    "empirical_es"
  )
  ls <- x[sub.vec]
  
  # Add ses and es suffixes to column names
  colnames(ls$empirical_es) <- paste0(colnames(ls$empirical_es),"_es")
  
  ## use this code if using ses rather than es
  # colnames(ls$ses) <- paste0(colnames(ls$ses),"_ses")  
  
  # Merge obs, es, and ses of each measure
  ls <- lapply(ls, function(y) {
    as.data.frame(y) %>% tibble::rownames_to_column("COMID")
    })
  purrr::reduce(ls,left_join,by=join_by(COMID)) 
})

# Export
saveRDS(final_list$alpha,file.path(export_dir,"native_alpha.rds"))
saveRDS(final_list$his_lcbd,file.path(export_dir,"native_lcbd.rds"))
saveRDS(final_list$`_d_lcbd`,file.path(export_dir,"delta_lcbd.rds"))


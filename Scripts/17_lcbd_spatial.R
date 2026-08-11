#-------------------------------------------------------------------------------
#
#  Beta diversity change spatial plots  
#
#-------------------------------------------------------------------------------

# Author: 

# Created: 04/30/2026

# Description:
# This script aggregates changes in LCBD values to HUC6 resolution to visualize 
# spatial trends in delta LCBD across the United States.

 #House keeping ----------------------------------------------------------------
rm(list = ls())

# Directories
beta_dir <- "Diversity Output Data"
pred_dir <- "Analysis_data"
export_dir <- "Results"

# Load packages
library(sf)
library(tidyr)
library(dplyr)

# Load data
lcbd_df <- readRDS(file.path(beta_dir,"delta_lcbd.rds"))
inv_df <-readRDS(file.path(pred_dir,"origin_invaded.rds"))
huc_df <- st_read(file.path(pred_dir,"huc6.shp"))


# Prepare data -----------------------------------------------------------------

# Bridge COMID to HUC6
lcbd_huc <- inv_df %>% 
  mutate(HUC6 = substr(HUC_12,1,6)) %>%   
  select(HUC6,COMID) %>% 
  right_join(lcbd_df %>% mutate(COMID = as.numeric(COMID))) %>% 
  select(!contains("ses")) %>% 


# Calculate the mean LBD changes in each huc6
  group_by(HUC6) %>%  
  summarise(
    n = n_distinct(COMID),
    across(-COMID, ~mean(.x))
  )
          
# Merge data into shape files
spatial_data <- huc_df %>% right_join(lcbd_huc)


# Prepare export ---------------------------------------------------------------

# Find min and max raw beta change values for plot labels
lcbd_min <-lcbd_huc %>% select(tax_Btotal,fun_Btotal,phy_Btotal) %>% max()
lcbd_max <-lcbd_huc %>% select(tax_Btotal,fun_Btotal,phy_Btotal) %>% min()

plot_data <- spatial_data %>% 
  mutate(
    plot_range = seq(lcbd_min,lcbd_max,length.out =nrow(spatial_data)),  # allows for even color scheme among plots
    plot_range_es = seq(-3.08,3.08,length.out =nrow(spatial_data))
  )  

# Export
st_write(plot_data,file.path(export_dir,"lcbd_map_data.shp"),delete_layer=T)

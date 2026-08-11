#-------------------------------------------------------------------------------
#
#  Beta diversity RDA predictor preparation
# 
#-------------------------------------------------------------------------------

# Author: 

# Created:04/28/2026

# Description: Compiles and formats community invadedness, propagule pressure, 
# habitat characteristic, habitat alteration, and biotic explanatory variables 
# for beta diversity change RDA and variance partitioning. Then community 
# invadedness summary stats are calculated.


# House keeping ----------------------------------------------------------------
rm(list = ls())

# Directories
com_dir <- "Diversity Input Data"
div_dir <- "Diversity Output Data"
pred_dir <- "Analysis_data"

# Load packages
library(dplyr)
library(tidyr)
# library(vegan)
library(StreamCatTools)

# Load data
inv_df <- readRDS(file.path(pred_dir,"origin_invaded.rds"))
alpha_df <-readRDS(file.path(div_dir,"native_alpha.rds"))
lcbd_df <- readRDS(file.path(div_dir,"native_lcbd.rds"))
fishing_df <- read.csv(file.path(pred_dir,"FreshwaterFishing_RecreationDemand.csv"))
hai_df<- read.csv(file.path(pred_dir,"predicted_alteration_region_model.csv"))
comid <- unique(inv_df$COMID)

# Origin-based invadedness -----------------------------------------------------

# Select only columns for RDA analysis
inv_for <- inv_df %>% 
  select(COMID,invProv,invReg,invExtr,nTaxA)


# Native community structure  -------------------------------------------------
nat_df <- alpha_df %>% 
  left_join(lcbd_df) %>% 
  select(!contains("ses")) %>% 
  rename(
    nFunA = fun_his_alpha,
    nPhyA = phy_his_alpha,
    nFunAES = fun_his_alpha_es,
    nPhyAES = phy_his_alpha_es,
    nTaxBtot = tax_Btotal,
    nTaxBrep = tax_Brepl,
    nTaxBric = tax_Brich,
    nFunBtot = fun_Btotal,
    nFunBrep = fun_Brepl,
    nFunBric = fun_Brich,
    nPhyBtot = phy_Btotal,
    nPhyBrep = phy_Brepl,
    nPhyBric = phy_Brich,
    nFunBtotES = fun_Btotal_es,
    nFunBrepES = fun_Brepl_es,
    nFunBricES = fun_Brich_es,
    nPhyBtotES = phy_Btotal_es,
    nPhyBrepES = phy_Brepl_es,
    nPhyBricES = phy_Brich_es
  ) %>% 
  mutate(
    COMID = as.numeric(COMID),
    "nTaxRepRic" = nTaxBrep/nTaxBric,
    "nFunRepRic" = nFunBrep/nFunBric,
    "nPhyRepRic" = nPhyBrep/nPhyBric
    )


# Propagule pressure  ---------------------------------------------------------
prop_df <- fishing_df %>% 
  
  # Change to 12 digit charater with leading zeros as needed
  mutate(
    HUC_12 = stringr::str_pad(HUC_12, width = 12, side = "left", pad = "0")
    ) %>% 
  
  # Assign fishing demand to COMIDs
  right_join(inv_df %>% distinct(COMID,HUC_12), by=join_by(HUC_12)) %>% 
  rename(FishDemand = FF_Demand) %>% 
  select(-HUC_12)


# Streamcat abiotic factors ----------------------------------------------------

# Metrics to download
sc_metrics <- c(
  "elev",
  "tmean9120",
  "precip9120",
  "bfi"
  )

# Dowload and format streamcat data
sc_df <-StreamCatTools::sc_get_data(
    comid =comid, 
    metric = sc_metrics,
    aoi="ws",
    showAreaSqKm = T) %>% 
  rename(
    COMID = comid,
    wsArea = wsareasqkm,
    wsElev = elevws,
    wsTemp = tmean9120ws,
    wsPrecip = precip9120ws,
    wsBFI = bfiws
  ) %>% 
  select(-catareasqkm,-catareasqkmrp100,-wsareasqkmrp100)


# Hydrological alteration ------------------------------------------------------
hai_for <- hai_df %>%

  # Format and subset to only include data in analysis
  select(COMID,pnHA_rank) %>%
  rename(HAI = pnHA_rank) %>%
  filter(COMID %in% comid)


# Compilation ------------------------------------------------------------------
predictor_df <- inv_for %>% 
  left_join(nat_df) %>% 
  left_join(prop_df) %>% 
  left_join(sc_df) %>% 
  left_join(hai_for) %>% 
  arrange(COMID) %>% 
  tidyr::drop_na() 

# CHeck for correlation
cor(predictor_df %>% select(-COMID))


# Export -----------------------------------------------------------------------
saveRDS(predictor_df,file.path(pred_dir,"rda_predictors.rds"))


# Invadedness Summary stats ----------------------------------------------------

# Mean invadedness 
invd_summary_data <-inv_df %>% 
  mutate(prop_nn = (prov + reg + extr)/tot_sp) # Create a columne for proportion invaded for all orgins
summary(invd_summary_data)

# standard deviations
sd(invd_summary_data$invProv)
sd(invd_summary_data$invReg)
sd(invd_summary_data$invExtr)


#  Proportion of sites with at least one nonnative species of each origin 
nrow(invd_summary_data[invd_summary_data$reg > 0,])/nrow(invd_summary_data)   # Regional
nrow(invd_summary_data[invd_summary_data$extr > 0,])/nrow(invd_summary_data)  # Extra realm
nrow(invd_summary_data[invd_summary_data$prov > 0,])/nrow(invd_summary_data)  # provincially
nrow(invd_summary_data[invd_summary_data$prov > 0 |                           # any nonnative
                         invd_summary_data$reg > 0 |
                         invd_summary_data$extr > 0 ,])/nrow(invd_summary_data)



# Nonnative species tables -----------------------------------------------------

# tables of species by nonnative origin, and number of sites they occupy as
# that type of nonnative species
nn_table_df <- com_df %>% 
  mutate(
    native_status = case_when(
    Native8 == T ~ "native",  # native to that huc 8
    Native8 == F & Native2 == T ~ "prov",  # prov species, native to region but not to huc 8
    Native2 ==F & NativeCon ==T ~ "reg",  # regional species, native continent but not regions
    NativeCon == F ~ "extr" # extra realm, not native to continent
  )
) 

# extra-realm species
sp_summary_exotic <- nn_table_df %>%   # summarize by cluster
  filter(native_status == "extr") %>% 
  group_by(Scientific_Name) %>% 
  summarise(n_sites = n()) %>% 
  mutate(prop_sites = n_sites/1023)

# Regional
sp_summary_reg <- nn_table_df %>%   # summarize by cluster
  filter(native_status == "reg") %>% 
  group_by(Scientific_Name) %>% 
  summarise(n_sites = n())%>% 
  mutate(prop_sites = n_sites/1023)

# Provincial
sp_summary_tran <- nn_table_df %>%   # summarize by cluster
  filter(native_status == "prov") %>% 
  group_by(Scientific_Name) %>% 
  summarise(n_sites = n())%>% 
  mutate(prop_sites = n_sites/1023)

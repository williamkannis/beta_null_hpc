#-------------------------------------------------------------------------------
#
#  Beta diversity change driver analysis
# 
#-------------------------------------------------------------------------------

# Author: 

# Created: 04/28/2026

# Description: Conducts redundancy analysis of the effects of species origin and
# PAB invasion drivers on changes in Local Contributions to Beta Diversity (LCBD).
# The importance of each variable is evaluated using variance partitioning. All
# analyses are repeated using raw LCBD values and null model standardized LCBD
# values. Plots and stats related to analyses are exported.


# House keeping ----------------------------------------------------------------
rm(list = ls())

# Directories
beta_dir <- "Diversity Output Data"
pred_dir <- "Analysis_data"
result_dir <- "Results"

# Load packages
library(tidyverse)  # maybe replace with ind pacakges
library(vegan)
# library(ggvegan)

# Load data
response_df <- readRDS(file.path(beta_dir,"delta_lcbd.rds"))
predictor_df <- readRDS(file.path(pred_dir,"rda_predictors.rds"))


# Prepare data -----------------------------------------------------------------

# Scale and center predictors
predictor_stand <- predictor_df %>% 
  tibble::column_to_rownames("COMID") %>% 
  decostand(method='standardize')

# Select only rows with predictor data
response_sub <- response_df %>% 
  select(!contains("ses")) %>% 
  filter(COMID %in% predictor_df$COMID) %>% 
  arrange(as.numeric(COMID)) %>% 
  tibble::column_to_rownames("COMID") %>% 
  rename(
    TaxBtot = tax_Btotal,
    TaxBrep = tax_Brepl,
    TaxBric = tax_Brich,
    FunBtot = fun_Btotal,
    FunBrep = fun_Brepl,
    FunBric = fun_Brich,
    PhyBtot = phy_Btotal,
    PhyBrep = phy_Brepl,
    PhyBric = phy_Brich,
    FunBtotES = fun_Btotal_es,
    FunBrepES = fun_Brepl_es,
    FunBricES = fun_Brich_es,
    PhyBtotES = phy_Btotal_es,
    PhyBrepES = phy_Brepl_es,
    PhyBricES = phy_Brich_es
  )
  


# Select predictors for observed and effect size rda
pred <- c("invProv","invReg","invExtr","nTaxA", "nFunAES", "nPhyAES","nTaxBtot",
           "FishDemand", "HAI", "wsElev", "wsArea","wsTemp","wsBFI")
pred_obs <- c(pred,"nFunBtot", "nPhyBtot")
pred_es <- c(pred,"nFunBtotES", "nPhyBtotES")

predictors_obs <- predictor_stand[,pred_obs]
predictors_es <- predictor_stand[,pred_es]

# Select responses for observed and effect size rda
res <- colnames(response_sub)
res_obs <- res[!grepl("ES",res)]
res_es <- c(
  res[grepl("ES",res) & !grepl("Tax",res)],
  res[!grepl("ES",res) & grepl("Tax",res)]
)

# Subset data based on raw or effect size
response_obs <- response_sub[,res_obs]
response_es <- response_sub[,res_es]

# Check to make sure rownames mathch
all(row.names(response_obs) == row.names(predictors_obs))
all(row.names(response_es) == row.names(predictors_es))

# Conduct initial rdas ----------------------------------------------------------
rda_init_obs <- rda(response_obs~., data = predictors_obs,scale=T)
rda_init_es <- rda(response_es~., data = predictors_es,scale=T)

# Check obs plots
plot(rda_init_obs, type = "n") 
points(rda_init_obs, pch=19, display = "sites") 
text(rda_init_obs, display = "species", col="blue") 
text(rda_init_obs, display = "bp", col="red") 

# Check ef plots
plot(rda_init_es, type = "n") 
points(rda_init_es, pch=19, display = "sites") 
text(rda_init_es, display = "species", col="blue") 
text(rda_init_es, display = "bp", col="red") 

# Variation explained by each axis
100 * rda_init_obs$CCA$eig / sum(rda_init_obs$CCA$eig)
100 * rda_init_es$CCA$eig / sum(rda_init_es$CCA$eig)

# Check vif
vif.cca(rda_init_obs)
vif.cca(rda_init_es)

#forward selection for obs
fs_obs <- ordistep(rda(response_obs~1, scale=T, data = predictors_obs), 
                scope = formula(rda_init_obs),  
                direction = "forward") 
fsR2_obs <- ordiR2step(rda(response_obs~1, scale=T, data = predictors_obs), 
                    scope = formula(rda_init_obs),  
                    direction = "forward")
# view results
fs_obs$anova
fsR2_obs$anova

#forward selection for es
fs_es <- ordistep(rda(response_es~1, scale=T, data = predictors_es), 
                   scope = formula(rda_init_es),  
                   direction = "forward") 
fsR2_es <- ordiR2step(rda(response_es~1, scale=T, data = predictors_es), 
                       scope = formula(rda_init_es),  
                       direction = "forward")
# View results
fs_es$anova
fsR2_es$anova

# Variables to retain for obs
fs_sig_obs <- sub("+ ","",row.names(fs_obs$anova),fixed = TRUE)
fsR2_sig_obs <- sub("+ ","",row.names(fsR2_obs$anova),fixed = TRUE)
sig_obs <- unique(c(fs_sig_obs,fsR2_sig_obs))
sig_obs <- sig_obs[sig_obs != "<All variables>"]

# Variables to retain for es
fs_sig_es <- sub("+ ","",row.names(fs_es$anova),fixed = TRUE)
fsR2_sig_es <- sub("+ ","",row.names(fsR2_es$anova),fixed = TRUE)
sig_es <- unique(c(fs_sig_es,fsR2_sig_es))
sig_es <- sig_es[sig_es != "<All variables>"]

# Subset precitors
pred_sig_obs <- predictors_obs[,sig_obs]
pred_sig_es <- predictors_es[,sig_es]

# predictors removed
setdiff(colnames(predictors_obs),sig_obs)
setdiff(colnames(predictors_es),sig_es)


# Conduct final rda -----------------------------------------------------------
rda_obs <- rda(response_obs~., data = pred_sig_obs,scale=T)
rda_es <- rda(response_es~., data = pred_sig_es,scale=T)

# Check obs plots
plot(rda_obs, type = "n") 
points(rda_obs, pch=19, display = "sites") 
text(rda_obs, display = "species", col="blue") 
text(rda_obs, display = "bp", col="red") 

# Check ef plots
plot(rda_es, type = "n") 
points(rda_es, pch=19, display = "sites") 
text(rda_es, display = "species", col="blue") 
text(rda_es, display = "bp", col="red") 

# Variation explained by each axis
100 * rda_obs$CCA$eig / sum(rda_obs$CCA$eig)
100 * rda_init_es$CCA$eig / sum(rda_init_es$CCA$eig)

# Importance tables
obs_important <-summary(rda_obs)[["cont"]][["importance"]]
es_important <-summary(rda_es)[["cont"]][["importance"]]

# axes loadings
sp_scor_obs <- as.data.frame(scores(rda_obs, display = "species"))
bp_scor_obs <- as.data.frame(scores(rda_obs, display = "bp"))
sp_scor_es <- as.data.frame(scores(rda_es, display = "species"))
bp_scor_es <- as.data.frame(scores(rda_es, display = "bp"))


# Significance of terms
anova(rda_obs, by="term") #not the best, results depend on the order of terms
anova(rda_obs, by="margin") #perform separate significance test for each marginal term in a model with all other terms
anova(rda_es, by="term") #not the best, results depend on the order of terms
anova(rda_es, by="margin") #perform separate significance test for each marginal term in a model with all other terms

# Check vif
max(vif.cca(rda_obs))
max(vif.cca(rda_init_es))

# Obs model fit
anova.cca(rda_obs, step = 1000)
obs_term_anova <- anova.cca(rda_obs, step = 1000, by = "term")
anova.cca(rda_obs, step = 1000, by = "axis")
RsquareAdj(rda_obs)

# Es model fit
anova.cca(rda_es, step = 1000)
es_term_anova <- anova.cca(rda_es, step = 1000, by = "term")
anova.cca(rda_es, step = 1000, by = "axis")
RsquareAdj(rda_es)

# Export RDA tables ------------------------------------------------------------

# Term significance (Table s8.1)
write.csv(obs_term_anova, file.path(result_dir,"obs_term_anova.csv"))
write.csv(es_term_anova,file.path(result_dir,"es_term_anova.csv"))

# Eigenvectors (Table s8.2)
write.csv(obs_important,file.path(result_dir,"obs_important.csv"))
write.csv(es_important,file.path(result_dir,"es_important.csv"))

# Loadings (Table s 8.3-4)
write.csv(bp_scor_obs,file.path(result_dir,"bp_scor_obs.csv"))
write.csv(sp_scor_obs,file.path(result_dir,"sp_scor_obs.csv"))
write.csv(bp_scor_es,file.path(result_dir,"bp_scor_es.csv"))
write.csv(sp_scor_es,file.path(result_dir,"sp_scor_es.csv"))


# RDA plots (Fig. 6)  -------------------------------------------------------------------

# Observed plots
png(
  filename = file.path(result_dir,"obs_rda_temp.png"),
  height = 1600,
  width = 2000)
par(mar=c(8,8,4,2)+1.2)
par(mgp=c(6,2,0))
plot(rda_obs, type = "n",
     xlab = "",
     ylab = "",
     xlim =c(-4,4),
     cex = 2.5,
     cex.lab = 2.5,
     cex.axis = 3,
     lwd=1)
abline(h = 0, v = 0, lwd = 2,lty=5)
points(rda_obs, pch=19, display = "species",cex = 4,col="darkgrey") 
text(rda_obs, display = "bp", col="red", cex = 2.5)
dev.off()

png(
  filename = file.path(result_dir,"obs_rda_blank.png"),
  height = 1600,
  width = 2000)
par(mar=c(8,8,4,2)+1.2)
par(mgp=c(6,2,0))
plot(rda_obs, type = "n",
     xlab = "",
     ylab = "",
     xlim =c(-4,4),
     cex = 2.5,
     cex.lab = 2.5,
     cex.axis = 3,
     lwd=1)
abline(h = 0, v = 0, lwd = 2,lty=5)
points(rda_obs, pch=19, display = "sites",cex = 4,col="darkgrey") #we want to show sites scores as points
dev.off()

# es plots
png(
  filename = file.path(result_dir,"es_rda_temp.png"),
  height = 1600,
  width = 2000)
par(mar=c(8,8,4,2)+1.2)
par(mgp=c(6,2,0))
plot(rda_es, type = "n",
     xlab = "",
     ylab = "",
     xlim =c(-2,3),
     cex = 2.5,
     cex.lab = 2.5,
     cex.axis = 3,
     lwd=1)
abline(h = 0, v = 0, lwd = 2,lty=5)
points(rda_es, pch=19, display = "species",cex = 4,col="darkgrey") 
text(rda_es, display = "bp", col="red", cex = 2.5)
dev.off()

png(
  filename = file.path(result_dir,"es_rda_blank.png"),
  height = 1600,
  width = 2000)
par(mar=c(8,8,4,2)+1.2)
par(mgp=c(6,2,0))
plot(rda_es, type = "n",
     xlab = "",
     ylab = "",
     xlim =c(-2,3),
     cex = 2.5,
     cex.lab = 2.5,
     cex.axis = 3,
     lwd=1)
abline(h = 0, v = 0, lwd = 2,lty=5)
points(rda_es, pch=19, display = "sites",cex = 4,col="darkgrey") #we want to show sites scores as points
dev.off()

  
# Variance partitioning  ---------------------------------------------------------

# Set vars
inv_vars <- c("invExtr","invProv","invReg")
div_vars <- c("nTaxBtot","nFunBtot","nPhyBtot","nTaxA","nFunA","nPhyA",
              "nFunBtotES","nPhyBtotES","nTaxAES","nFunAES","nPhyAES")
prop_vars <- c("FishDemand")
ab_vars <- c("HAI", "wsElev", "wsArea","wsTemp","wsBFI")

# Subset observed data
obs_inv <- pred_sig_obs %>% select(contains(inv_vars))
obs_div <- pred_sig_obs %>% select(contains(div_vars))
obs_prop <- pred_sig_obs %>% select(contains(prop_vars))
obs_ab <- pred_sig_obs %>% select(contains(ab_vars))

# Var partitioning of observed data
obs_varpart <- varpart(decostand(response_obs, method = "standardize"),obs_inv,obs_div,obs_ab,obs_prop)
plot(obs_varpart)
obs_varpart$part

# Subset effect size data
es_inv <- pred_sig_es %>% select(contains(inv_vars))
es_div <- pred_sig_es %>% select(contains(div_vars))
es_prop <- pred_sig_es %>% select(contains(prop_vars))
es_ab <- pred_sig_es %>% select(contains(ab_vars))

# Var partitioning of effect size data
es_varpart <- varpart(decostand(response_es, method = "standardize"),es_inv,es_div,es_ab,es_prop)
plot(es_varpart)
es_varpart$part


# Subset data for significance testing
obs_pab <- pred_sig_obs %>% select(contains(c(prop_vars,ab_vars,div_vars)))
obs_ipa <- pred_sig_obs %>% select(contains(c(prop_vars,ab_vars,inv_vars)))
obs_ipb <- pred_sig_obs %>% select(contains(c(prop_vars,inv_vars,div_vars)))
obs_iab <- pred_sig_obs %>% select(contains(c(inv_vars,ab_vars,div_vars)))

es_pab <- pred_sig_es %>% select(contains(c(prop_vars,ab_vars,div_vars)))
es_ipa <- pred_sig_es %>% select(contains(c(prop_vars,ab_vars,inv_vars)))
es_ipb <- pred_sig_es %>% select(contains(c(prop_vars,inv_vars,div_vars)))
es_iab <- pred_sig_es %>% select(contains(c(inv_vars,ab_vars,div_vars)))

# effect size Significance
anova.cca(rda(decostand(response_obs, method = "standardize"),obs_inv,obs_pab))
anova.cca(rda(decostand(response_obs, method = "standardize"),obs_div,obs_ipa))
anova.cca(rda(decostand(response_obs, method = "standardize"),obs_ab,obs_ipb))
anova.cca(rda(decostand(response_obs, method = "standardize"),obs_prop,obs_iab))

# effect size Significance
anova.cca(rda(decostand(response_es, method = "standardize"),es_inv,es_pab))
anova.cca(rda(decostand(response_es, method = "standardize"),es_div,es_ipa))
anova.cca(rda(decostand(response_es, method = "standardize"),es_ab,es_ipb))
anova.cca(rda(decostand(response_es, method = "standardize"),es_prop,es_iab))


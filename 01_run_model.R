library(patchwork)
library(splm)
library(plm)
library(tidyverse)
library(sf)
library(spatialreg)
library(brazilmaps)

# Load data formatted for panel analysis
panel_data <- read_csv("data/panel_amazon_lag.csv")

# Drop municiaplities with no deforestation in every time period
nodeforest <- 
    panel_data %>% 
    filter(avg_yearly_deforestation == 0 | lagged_deforestation == 0) %>% 
    count(ibge7_code) %>% 
    filter(n>0) %>% 
    select(ibge7_code)

panel_data <- 
    panel_data %>% 
    filter(!(ibge7_code %in% nodeforest$ibge7_code))

# Format non-numeric data as factors for model
panel_data  <- panel_data %>%
    mutate_at(1, as.character) %>%
    mutate(GE_consortiums = case_when(panel_data$GE_consortiums != 0 ~ 1,
                                      TRUE ~ 0)) %>%
    mutate(ROL_U_improvement = case_when(panel_data$ROL_U_improvement  != 0 ~ 1,
                                         TRUE ~ 0)) %>%
    mutate(ROL_U_neighborhood = case_when(panel_data$ROL_U_neighborhood  != 0 ~ 1,
                                          TRUE ~ 0)) %>%
    mutate(ROL_zoning = case_when(panel_data$ROL_zoning != 0 ~ 1,
                                  TRUE ~ 0)) %>%
    mutate(GE_masterplan = case_when(panel_data$GE_masterplan != 0 ~ 1,
                                     TRUE ~ 0)) %>%
    mutate_at(c("EG_envir.council", "EG_envir.fund", "EG_envir.agency", 
                "GE_consortiums", "RQ_enterprise_incentive_existence", 
                "RQ_enterprise_restriction_existence", "VA_webpage", "VA_mayor_woman",
                "state", "GE_masterplan","ROL_U_improvement",
                "ROL_U_neighborhood", "ROL_zoning", "GE_masterplan"), as_factor) 

# Drop one municipality that doesn't have complete records on RQ enterprice restrictions
panel_data <- panel_data[panel_data$ibge7_code != "1400050", ]
panel_data <- arrange(panel_data, ibge7_code)

# Get spatial polygons for municipalities
municipalities <-
    get_brmap(geo = "City") %>%
    filter(State %in% c(11,12,13,14,15,16,17,21,51)) %>%
    rename(ibge7_code = "City") %>%
    filter(ibge7_code %in% panel_data$ibge7_code) 

states <- 
    get_brmap(geo="State") %>%
    filter(State %in% c(11,12,13,14,15,16,17,21,51)) 

municipalities$ibge7_code <- as.character(municipalities$ibge7_code)

# Build spatial neighbords list and spatial weights
muni_sp <- municipalities %>% as("Spatial")
nb <- poly2nb(muni_sp, queen=TRUE)
lw <- nb2listw(nb, style="W", zero.policy=NULL)

## Prepare subsets of data for model evaluation
panel_data_short <- panel_data %>% filter(snapshot != "TS4")
panel_data_tmp <- panel_data %>% filter(snapshot != "TS1")
panel_data_tmp_vars <- panel_data_short %>% select(EG_InterMunicCop:pa)
panel_data_splag <- tibble(select(panel_data_tmp, 
                                  ibge7_code:lagged_rel_forest_area), 
                                  panel_data_tmp_vars)

# Model Tests -----------
avg_fixed_fm <- log(avg_yearly_deforestation) ~ 
    log(lagged_deforestation) + crop_dens + cattle_dens + popden + GDP_reais +
    EG_envir.council + EG_envir.agency + EG_envir.fund + EG_envir.person + 
    GE_employees + GE_consortiums + GE_masterplan + 
    ROL_div_land + ROL_U_improvement + ROL_U_neighborhood + ROL_zoning + 
    RQ_ag_comp + RQ_nonag_comp + RQ_ag_personnel + RQ_non_ag_personnel + 
    RQ_enterprise_incentive_existence + RQ_enterprise_restriction_existence + 
    VA_n_candidates + VA_n_compan_communic + VA_prop_votes + 
    VA_webpage + VA_mayor_woman + factor(snapshot)

# Alternate formula 1 - controls only
fm_nogov <- log(avg_yearly_deforestation) ~ 
    log(lagged_deforestation) + crop_dens + cattle_dens + popden + GDP_reais + factor(snapshot)

# Alternate formula 2 - Sig. only
fm_sig <- log(avg_yearly_deforestation) ~ 
    log(lagged_deforestation) + crop_dens + cattle_dens + popden  +
    + GDP_reais + factor(snapshot) + EG_envir.agency + EG_envir.fund +
    RQ_ag_personnel + RQ_non_ag_personnel + VA_mayor_woman 

# Alternate formula 3 - EG/RQ
fm_eg_rq <- log(avg_yearly_deforestation) ~ 
    log(lagged_deforestation) + crop_dens + cattle_dens + popden +
    + GDP_reais +  EG_envir.council + EG_envir.agency + EG_envir.fund + 
    EG_envir.person + factor(snapshot) + RQ_ag_comp + RQ_nonag_comp + 
    RQ_ag_personnel + RQ_non_ag_personnel + 
    RQ_enterprise_incentive_existence + RQ_enterprise_restriction_existence


# Model S1 - Preferred model (lagged, spatial)
mod_s1 <- spml(formula = avg_fixed_fm, 
                             data = panel_data_splag, 
                             index=c("ibge7_code", "snapshot"), 
                             listw = lw, model = "within", 
                             lag = F,
                             spatial.error = "kkp", 
                             zero.policy=TRUE)

# Model S2 - Controls only
mod_s2 <- spml(formula = fm_nogov, 
                         data = panel_data_splag, 
                         index=c("ibge7_code", "snapshot"), 
                         listw = lw, 
                         model = "within", 
                         lag = F,
                         spatial.error = "kkp", 
                         zero.policy=TRUE)


# Model S3 - Significant variables only
mod_s3 <- spml(formula = fm_sig, 
                         data = panel_data_splag, 
                         index=c("ibge7_code", "snapshot"), 
                         listw = lw, model = "within", 
                         lag = F,
                         spatial.error = "kkp", 
                         zero.policy=TRUE)

# Model S4 - EG and RQ only
mod_s4 <- spml(formula = fm_eg_rq, 
               data = panel_data_splag, 
               index=c("ibge7_code", "snapshot"), 
               listw = lw, model = "within", 
               lag = F,
               spatial.error = "kkp", 
               zero.policy=TRUE)

# Model S5 - Unlagged
mod_s5 <- spml(formula = avg_fixed_fm, 
                         data = panel_data_short, 
                         listw = nb2listw(nb, style="W", zero.policy=TRUE),
                         model = "within", 
                         lag = F,
                         spatial.error = "kkp", 
                         zero.policy=TRUE)



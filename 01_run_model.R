library(patchwork)
library(scales)
library(splm)
library(plm)
library(tidyverse)
library(sf)
library(spatialreg)
library(brazilmaps)
library(gridExtra)
library(grid)
library(ggfortify)


##
theme_Publication <- function(base_size=13) {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.2)),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =1),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                #legend.direction = "horizontal",
                legend.key.size= unit(0.5, "cm"),
                #legend.margin = unit(0, "cm"),
                legend.title = element_blank(),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
}

# load data formatted for panel analysis
panel_data <- read_csv("panel_amazon_lag.csv")

# drop municiaplities with no deforestation in every time period
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
panel_data_splag <- tibble(select(panel_data_tmp, ibge7_code:lagged_rel_forest_area), panel_data_tmp_vars)

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

# Fixed effects
avg_fixed_mod <- plm(formula = avg_fixed_fm, data=panel_data, 
                     model="within")

# Lagged fixed effects
avg_fixed_mod_lag <- plm(formula = avg_fixed_fm, data=panel_data_splag, 
                         model="within", 
                         index=c("ibge7_code", "snapshot"))

# Spatial fixed effects
avg_fixed_mod_sp <- spml(formula = avg_fixed_fm, data = panel_data_short, 
                         listw = nb2listw(nb, style="W", zero.policy=TRUE), 
                         model = "within", lag = F,
                         spatial.error = "kkp", zero.policy=TRUE)

# Lagged spatial fixed effects
avg_fixed_mod_lag_sp <- spml(formula = avg_fixed_fm, data = panel_data_splag, 
                             index=c("ibge7_code", "snapshot"), 
                             listw = lw, model = "within", lag = F,
                             spatial.error = "kkp", zero.policy=TRUE)

# Alternate formula for robustness checks
fm_nogov <- log(avg_yearly_deforestation) ~ 
    log(lagged_deforestation) + crop_dens + cattle_dens + popden + GDP_reais + factor(snapshot)

# Fixed effects no governance variables
avg_fixed_mod_ng <- plm(formula = fm_nogov, data=panel_data, 
                        model="within")

# Lagged fixed effects no governance variables
avg_fixed_mod_lag_ng <- plm(formula = fm_nogov, data=panel_data_splag, 
                            model="within", 
                            index=c("ibge7_code", "snapshot"))

# Spatial fixed effects no governance variables
avg_fixed_mod_sp_ng <- spml(formula = fm_nogov, data = panel_data_short, 
                            listw = nb2listw(nb, style="W", zero.policy=TRUE),
                            model = "within", lag = F,
                            spatial.error = "kkp", zero.policy=TRUE)

# Spatial lagged fixed effects no governance variables

avg_fixed_mod_lag_sp_ng <- spml(formula = fm_nogov, data = panel_data_splag, 
                                index=c("ibge7_code", "snapshot"), listw = lw, model = "within", lag = F,
                                spatial.error = "kkp", zero.policy=TRUE)


# Model Tests -----------
snap_labels = c(TS1 = "2005-2008",
                TS2 = "2009-2012",
                TS3 = "2013-2016",
                TS4 = "2017-2018")

# Deforestation plot
left_join(municipalities, panel_data) %>%
    ggplot() + 
    geom_sf(aes(fill=avg_yearly_deforestation-lagged_deforestation), lwd=0.0) + 
    geom_sf(data=states, fill=NA, lwd=0.5) + 
    scale_fill_gradient2(limits = c(-250,150), oob=squish, low=muted("blue"), mid="lightgrey",
                         high=muted("red"),
                         name="Avg. Yearly Deforestation Change (Sq. Km)",
                         labels=c("<-225", "-150", "-75", "0", "75", ">150"), 
                         breaks=c(-225, -150, -75, 0, 75, 150),
                         guide= guide_colourbar(title.position="top",barwidth=10,
                                                title.hjust=1)) +
    #ggtitle("Change in deforestation from previous period") + 
    facet_wrap(~ snapshot, labeller = labeller(snapshot=snap_labels)) +
    theme_minimal() +  
    theme(plot.title = element_text(size = 15),
          legend.position = "bottom",
          axis.text = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank())  

#ggsave("map1.png", width = 8, height = 6)

######## Residuals
gen_resid_spat_fe <- function(model) {
    # Formats residuals for spatial object
    df <- data.frame(x = names(model$residuals))
    resid <- df %>% 
        separate(x, c("ibge7_code", "snapshot")) %>%
        mutate(resid = model$residuals)
    return(resid)
}

# Get residuals from model objects for plotting
plm_fit <- panel_data_short %>% 
    dplyr::select(ibge7_code, snapshot, avg_yearly_deforestation) %>%
    mutate(fit_avg_fixed = 
               avg_fixed_mod$model$`log(avg_yearly_deforestation)` - 
               avg_fixed_mod$residuals,
           res = avg_yearly_deforestation-exp(fit_avg_fixed),
           res_plot = avg_fixed_mod$residuals)

splm_nolag_resid <- gen_resid_spat_fe(avg_fixed_mod_sp)
splm_fit <- panel_data_short %>% 
    arrange(snapshot, ibge7_code) %>% 
    dplyr::select(ibge7_code, snapshot, avg_yearly_deforestation) %>%
    mutate(fit_avg_fixed_sp = log(avg_yearly_deforestation) - avg_fixed_mod_sp$residuals,
           res = avg_fixed_mod_sp$residuals) %>%
    left_join(splm_nolag_resid)

plm_fit_lag <- panel_data_splag  %>% 
    dplyr::select(ibge7_code, snapshot, avg_yearly_deforestation) %>%
    mutate(fit_avg_fixed = 
               avg_fixed_mod_lag$model$`log(avg_yearly_deforestation)` - 
               avg_fixed_mod_lag$residuals,
           res = avg_yearly_deforestation-exp(fit_avg_fixed),
           res_plot = avg_fixed_mod_lag$residuals)

splm_lag_resid <- gen_resid_spat_fe(avg_fixed_mod_lag_sp)
splm_fit_lag <- panel_data_splag  %>% 
    arrange(snapshot, ibge7_code) %>% 
    dplyr::select(ibge7_code, snapshot, avg_yearly_deforestation) %>%
    mutate(fit_avg_fixed_sp = log(avg_yearly_deforestation) - 
               avg_fixed_mod_lag_sp$residuals,
           res = avg_fixed_mod_lag_sp$residuals) %>%
    left_join(splm_lag_resid)

# Residual plots
r1 <- ggplot(plm_fit)+
    geom_point(aes(fit_avg_fixed, res_plot, color=snapshot), alpha=0.5) +
    scale_color_manual(values=c("#440154FF", "#31688EFF", "#35B779FF"),
                       labels= c("2005-2008", "2009-2012", "2013-2016")) +
    geom_hline(yintercept=0, linetype=2, size=1, colour='darkgrey') + 
    ggtitle("Fixed Effects, No Lag") +
    theme(plot.title = element_text(size = 10)) +
    xlab("Fitted Value") +
    ylab("Residual") +
    theme_Publication()

r2 <- ggplot(splm_fit)+
    geom_point(aes(fit_avg_fixed_sp, resid, color=snapshot), alpha=0.5) +
    scale_color_manual(values=c("#440154FF", "#31688EFF", "#35B779FF"),
                       labels= c("2005-2008", "2009-2012", "2013-2016")) +
    geom_hline(yintercept=0, linetype=2, size=1, colour='darkgrey') + 
    ggtitle("(b) Residuals, Non-Lagged Model") +
    theme(plot.title = element_text(size = 10)) +
    xlab("Fitted Value (Log Scale)") +
    ylab("Residual") +
    theme_Publication()

f2 <- ggplot(splm_fit) +
    geom_point(aes(avg_yearly_deforestation, exp(fit_avg_fixed_sp), color=snapshot), alpha=0.5) + 
    scale_color_manual(values=c("#440154FF", "#31688EFF", "#35B779FF"),
                       labels= c("2005-2008", "2009-2012", "2013-2016")) +
    ggtitle("(a) Model Fit, Non-Lagged Model") +
    xlim(0,350) + ylim(0,400) + 
    xlab("Avg. Yearly Deforestation (sqkm/yr)") +
    ylab("Fitted Values") + 
    theme(plot.title = element_text(size = 10)) +
    theme_Publication()

r3 <- ggplot(plm_fit_lag)+
    geom_point(aes(fit_avg_fixed, res_plot, color=snapshot), alpha=0.5) +
    scale_color_manual(values=c("#31688EFF", "#35B779FF", "#FDE725FF"),
                       labels= c("2009-2012", "2013-2016", "2017-2018")) +
    geom_hline(yintercept=0, linetype=2, size=1, colour='darkgrey') + 
    ggtitle("Fixed Effects, Lag") +
    theme(plot.title = element_text(size = 10)) +
    xlab("Fitted Value") +
    ylab("Residual") +
    theme_Publication()

r4 <- ggplot(splm_fit_lag)+
    geom_point(aes(fit_avg_fixed_sp, resid, color=snapshot), alpha=0.5) +
    geom_hline(yintercept=0, linetype=2, size=1, colour='darkgrey') + 
    scale_color_manual(values=c("#31688EFF", "#35B779FF", "#FDE725FF"),
                       labels= c("2009-2012", "2013-2016", "2017-2018")) +
    ggtitle("Residuals, Lagged Model") +
    theme(plot.title = element_text(size = 10)) +
    xlab("Fitted Value (Log Scale)") +
    ylab("Residual") +
    theme_Publication()

f4 <- ggplot(splm_fit_lag) +
    geom_point(aes(avg_yearly_deforestation, exp(fit_avg_fixed_sp), color=snapshot), alpha=0.5) + 
    scale_color_manual(values=c("#31688EFF", "#35B779FF", "#FDE725FF"),
                       labels= c("2009-2012", "2013-2016", "2017-2018")) +
    ggtitle("Model Fit, Lagged Model") +
    xlim(0,350) + ylim(0,400) + 
    xlab("Avg. Yearly Deforestation (sqkm/yr)") +
    ylab("Fitted Values") + 
    theme(plot.title = element_text(size = 10)) +
    theme_Publication()

# Spatial Residuals --------------
snap_labels_nolag = c(TS1 = "2005-2008",
                TS2 = "2009-2012",
                TS3 = "2013-2016")

snap_labels_lag = c(TS2 = "2009-2012",
                    TS3 = "2013-2016",
                    TS4 = "2017-2018")

sp1 <- ggplot(left_join(municipalities, plm_fit)) + 
    geom_sf(aes(fill=res_plot), lwd=0.0) + 
    scale_fill_gradient2() +
    facet_wrap(~ snapshot, labeller = labeller(snapshot=snap_labels_nolag), drop=T) +
    ggtitle("Fixed Effects, Unlagged") +
    theme(plot.title = element_text(size = 10)) +
    theme_Publication()

sp2 <- ggplot(left_join(municipalities, plm_fit_lag)) + 
    geom_sf(aes(fill=res_plot), lwd=0.0) + 
    scale_fill_gradient2() +
    facet_wrap(~ snapshot, labeller = labeller(snapshot=snap_labels_lag), drop=T) +
    ggtitle("Fixed Effects, Lagged") +
    theme(plot.title = element_text(size = 10)) +
    theme_Publication()

sp3 <- ggplot(left_join(municipalities, splm_fit)) + 
    geom_sf(aes(fill=resid), lwd=0.0) + 
    scale_fill_gradient2() +
    facet_wrap(~ snapshot, labeller = labeller(snapshot=snap_labels_nolag), drop=T) +
    ggtitle("Spatial Fixed Effects, Unlagged") +
    theme(plot.title = element_text(size = 10)) +
    theme_Publication()
    
sp4 <- ggplot(left_join(municipalities, splm_fit_lag)) + 
    geom_sf(aes(fill=resid), lwd=0.0) + 
    scale_fill_gradient2() +
    facet_wrap(~ snapshot, labeller = labeller(snapshot=snap_labels_lag), drop=T) +
    ggtitle("Spatial Fixed Effects, Lagged") +
    theme(plot.title = element_text(size = 10)) +
    theme_Publication()


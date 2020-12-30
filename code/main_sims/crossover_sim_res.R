library(knitr)
library(kableExtra)
library(here)
library(tidyverse)
library(rsimsum)
setwd(here("jon_stuff/cross_sims_2yr"))

sim_settings <- 
    expand.grid(design = c("events", "time", "parallel"),
                ve = c("waning", "constant"),
                yr2 = c("same", "half"))

pars <- expand.grid(design = c("events", "time", "parallel"),
                    ve = c("waning", "constant"),
                    param = c("intcpt", "slope"),
                    yr2 = c("same", "half"),
                    model = c("const", "linear", "pspline"),
                    simulation = 1:1e4)

pars_cross <- 
    expand.grid(design = c("events", "time"),
                ve = c("waning", "constant"),
                param = c("intcpt", "slope"),
                yr2 = c("same", "half"),
                model = c("const", "linear", "pspline"),
                simulation = 1:1e4)

ests <- expand.grid(design = c("events", "time", "parallel"),
                    ve = c("waning", "constant"),
                    yr2 = c("same", "half"),
                    time = c(0.5, 1, 1.5, 2),
                    model = c("const", "linear", "pspline"),
                    simulation = 1:1e4)

par_list <- vector("list", nrow(sim_settings))
cross_pars_list <- vector("list", 8)
ve_list <- vector("list", nrow(sim_settings))
decay_list <- vector("list", nrow(sim_settings))
t_cross_list <- vector("list", nrow(sim_settings))
n_cross_list <- vector("list", nrow(sim_settings)) 

j = 0

for(s in seq_len(nrow(sim_settings))) {
    
    filename =
        paste0("cross_sim_", 
               sim_settings$ve[s], "_",
               sim_settings$design[s], "_",
               sim_settings$yr2[s], ".Rds")
    
    sim_res = readRDS(filename)
    
    # get n_cross and t_cross
    n_cross_list[[s]] <- sapply(sim_res, "[[", "n_cross")
    t_cross_list[[s]] <- sapply(sim_res, "[[", "t_cross")
    
    # subset to only the times and models of interest
    ve_sub = 
        do.call(rbind, lapply(sim_res, 
               function(x)
                   subset(x$VE_ests,
                          time %in% c(0.5, 1, 1.5, 2) & 
                          model %in% c("const", "linear", "pspline"))))
    
    decay_sub = 
        do.call(rbind, 
                lapply(sim_res, 
                       function(x) 
                           subset(x$VE_decays,
                                  time %in% c(0.5, 1, 1.5, 2) & 
                                      model %in% 
                                      c("const", "linear", "pspline"))))
    
    par_sub = 
        do.call(rbind,
                lapply(sim_res,
                       function(x)
                           subset(x$par_ests,
                                  model %in% c("const", "linear", "pspline"))))
    par_sub$param[par_sub$model != "const"] = c("intcpt", "slope")
    par_sub$param[par_sub$model == "const"] = "intcpt"
    par_sub$simulation = as.numeric(par_sub$simulation)
    par_sub$truth = as.numeric(par_sub$truth)
    par_sub$est = as.numeric(par_sub$est)
    par_sub$var = as.numeric(par_sub$var)
    
    inds = which(ests$design == sim_settings$design[s] & 
                     ests$ve == sim_settings$ve[s] & 
                     ests$yr2 == sim_settings$yr2[s])
    
    par_inds = which(pars$design == sim_settings$design[s] & 
                         pars$ve == sim_settings$ve[s] & 
                         pars$yr2 == sim_settings$yr2[s])
    
    # merge
    ve_list[[s]] <- 
        left_join(ests[inds,], ve_sub,
                  by = c("time", "model", "simulation"))
    
    decay_list[[s]] <- 
        left_join(ests[inds,], decay_sub,
                  by = c("time", "model", "simulation"))
    
    par_list[[s]] <-
        left_join(pars[par_inds,], par_sub,
                  by = c("param", "model", "simulation"))
    
    
    if(sim_settings$design[s] != "parallel") {
        
        par_cross_sub = 
            do.call(rbind,
                    lapply(sim_res, "[[", "par_ests_cross"))
        
        par_cross_sub$param[par_cross_sub$model != "const"] = c("intcpt", "slope")
        par_cross_sub$param[par_cross_sub$model == "const"] = "intcpt"
        
        par_cross_sub$simulation = as.numeric(par_cross_sub$simulation)
        par_cross_sub$truth = as.numeric(par_cross_sub$truth)
        par_cross_sub$est = as.numeric(par_cross_sub$est)
        par_cross_sub$var = as.numeric(par_cross_sub$var)
        
        par_cross_inds = 
            which(as.character(pars_cross$design) == sim_settings$design[s] & 
                      as.character(pars_cross$ve) == sim_settings$ve[s] &
                      as.character(pars_cross$yr2) == sim_settings$yr2[s])
        
        cross_pars_list[[s]] <-
            left_join(pars_cross[par_cross_inds,], par_cross_sub,
                      by = c("param", "model", "simulation"))
    }
}

# combine the results
res_ve = do.call(rbind, ve_list)
res_decay = do.call(rbind, decay_list)
res_pars = do.call(rbind, par_list)
res_pars_cross = do.call(rbind, cross_pars_list)

# summarize crossover times and event counts
cross_stats = 
    data.frame(
        sim_settings,
        n_cross_mean = sapply(n_cross_list, mean, na.rm = T),
        n_cross_sd   = sapply(n_cross_list, sd, na.rm = T),
        t_cross_mean = sapply(t_cross_list, mean, na.rm = T),
        t_cross_mean = sapply(t_cross_list, sd, na.rm = T)
    )

# summarize
ve_sum = 
    res_ve %>% 
    group_by(design, ve, yr2, time, model) %>% 
    summarize(emp_var = var(est),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

decay_sum = 
    res_decay %>% 
    group_by(design, ve, yr2, time, model) %>% 
    summarize(emp_var = var(est),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

par_sum = 
    res_pars %>% 
    group_by(design, ve, yr2, model, param) %>% 
    summarize(emp_var = var(est), 
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var)))) 

par_cross_sum = 
    res_pars_cross %>% 
    group_by(design, ve, yr2, model, param) %>% 
    summarize(emp_var = var(est), 
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

par_cross_sum$design = ifelse(par_cross_sum$design == "events", 
                                "Stop at 150 cases", "Stop at 1 year")

par_sum = 
    rbind(par_sum, par_cross_sum) %>% 
    arrange(ve, design, model, param)

par_sum = par_sum[complete.cases(par_sum),]

# subset to get the true data generating mechanism and the pspline
ve_sum = 
    ve_sum %>% 
    filter(model == "pspline" | model == "linear") %>% 
    arrange(ve, yr2, design, model, time) %>% 
    ungroup()

decay_sum = 
    decay_sum %>% 
    filter(model == "pspline" | model == "linear") %>% 
    arrange(ve, yr2, design, model, time) %>% 
    ungroup()

ve_sum = 
    ve_sum[,c("ve", "yr2", "design", "model", "time", "emp_var", "covg")]

decay_sum = 
    decay_sum[,c("ve", "yr2", "design", "model", "time", "emp_var", "covg")]

res_sum = full_join(ve_sum, decay_sum, c("ve", "yr2", "design", "model", "time"))

# recode some of the variables
res_sum = 
    res_sum %>% 
    mutate(design = 
               case_when(design == "events" ~ "Cross at 150 cases",
                         design == "time" ~ "Cross at 1 year",
                         TRUE ~ "Parallel trial"),
           model = 
               case_when(model == "pspline" ~ "P-spline",
                         model == "linear" ~ "log-linear",
                         model == "const" ~ "constant"))

par_sum = 
    par_sum %>% 
    mutate(design = 
               case_when(design == "events" ~ "Cross at 150 cases",
                         design == "time" ~ "Cross at 1 year",
                         design == "parallel" ~ "Parallel trial",
                         TRUE ~ design),
           model = 
               case_when(model == "pspline" ~ "P-spline",
                         model == "linear" ~ "log-linear",
                         model == "const" ~ "Constant VE"),
           param = 
               case_when(param == "intcpt" ~ "Intercept",
                         TRUE ~ "Linear trend"))

res_sum[,6:9] = round(res_sum[,6:9], digits = 3)
par_sum[,6:7] = round(par_sum[,6:7], digits = 3)

# pivot par_sum wider
par_sum = 
    par_sum %>% 
    pivot_wider(names_from = c("param"), values_from = c("emp_var", "covg"))
par_sum = par_sum[,c(1:4,5,7,6,8)]

# generate tables
tab_half <- 
    res_sum %>% 
    filter(yr2 == "half") %>%
    arrange(desc(ve)) %>% 
    select(3:9)

tab_half$design[-seq(1,48,by=8)] = " "
tab_half$model[-seq(1,48,by=4)] = " "

tab_half %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 3, "$VE(t)$" = 3, "$\\Delta VE(t)$" = 3)) %>%
    pack_rows("Vaccine efficacy constant at 75%", 1, 24) %>% 
    pack_rows("Vaccine efficacy wanes from 85% to 35% over 1.5 years", 25, 48)

tab_same <- 
    res_sum %>% 
    filter(yr2 != "half") %>%
    arrange(desc(ve)) %>% 
    select(3:9)

tab_same$design[-seq(1,48,by=8)] = " "
tab_same$model[-seq(1,48,by=4)] = " "

tab_same %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 3, "$VE(t)$" = 3, "$\\Delta VE(t)$" = 3)) %>%
    pack_rows("Vaccine efficacy constant at 75%", 1, 24) %>% 
    pack_rows("Vaccine efficacy wanes from 85% to 35% over 1.5 years", 25, 48)

# par table
par_sum$design = factor(par_sum$design, 
                        levels = c("Stop at 150 cases", "Stop at 1 year",
                                   "Cross at 150 cases", "Cross at 1 year",
                                   "Parallel trial"))

partab_half <- 
    par_sum %>% 
    filter(yr2 == "half") %>%
    arrange(desc(ve), design, model)

partab_half$design = as.character(partab_half$design)
partab_half$model = as.character(partab_half$model)

partab_half$design[-seq(1,30,by=3)] = " "
partab_half = subset(partab_half, !(ve == "waning" & model == "Constant VE"))

partab_half[,c(1,4:8)] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 1, 
                       "log-linear model" = 2, "P-spline model" = 2)) %>%
    pack_rows("Vaccine efficacy constant at 75%", 1, 15) %>% 
    pack_rows("Vaccine efficacy wanes from 85% to 35% over 1.5 years", 16, 25)

partab_same <- 
    par_sum %>% 
    filter(yr2 == "same") %>%
    arrange(desc(ve), design, model)

partab_same$design = as.character(partab_same$design)
partab_same$model = as.character(partab_same$model)

partab_same = subset(partab_same, !(ve == "waning" & model == "Constant VE"))
partab_same$design[c(seq(2,15,by=3), seq(3,15,by=3), seq(17,25,by=2))] = " "

partab_same[,c(1,4:8)] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 1, 
                       "log-linear model" = 2, "P-spline model" = 2)) %>%
    pack_rows("Vaccine efficacy constant at 75%", 1, 15) %>% 
    pack_rows("Vaccine efficacy wanes from 85% to 35% over 1.5 years", 16, 25)

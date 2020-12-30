library(knitr)
library(kableExtra)
library(here)
library(tidyverse)
library(rsimsum)
setwd(here("jon_stuff/cross_dur_2yr"))

sim_settings <- 
    expand.grid(ve = c("waning", "constant"),
                crossdur = c(2, 8))

pars <- expand.grid(ve = c("waning", "constant"),
                    param = c("intcpt", "slope"),
                    crossdur = c(2,8),
                    model = c("linear", "pspline"),
                    simulation = 1:1e4)

ests <- expand.grid(ve = c("waning", "constant"),
                    crossdur = c(2,8),
                    time = c(0, 0.5, 1, 1.5, 2),
                    model = c("const", "linear", "pspline"),
                    simulation = 1:1e4)

par_list <- vector("list", nrow(sim_settings))
ve_list <- vector("list", nrow(sim_settings))
decay_list <- vector("list", nrow(sim_settings))
t_cross_list <- vector("list", nrow(sim_settings))
n_cross_list <- vector("list", nrow(sim_settings)) 

for(s in seq_len(nrow(sim_settings))) {
    
    filename =
        paste0("cross_sim_", 
               sim_settings$ve[s], "_",
               sim_settings$crossdur[s], ".Rds")
    
    sim_res = readRDS(filename)
    
    # get n_cross and t_cross
    n_cross_list[[s]] <- sapply(sim_res, "[[", "n_cross")
    t_cross_list[[s]] <- sapply(sim_res, "[[", "t_cross")
    
    # subset to only the times and models of interest
    ve_sub = 
        do.call(rbind, lapply(sim_res, 
               function(x)
                   subset(x$VE_ests,
                          time %in% c(0, 0.5, 1, 1.5, 2) & 
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
                                  model %in% c("linear", "pspline"))))
    par_sub$param = c("intcpt", "slope")
    par_sub$simulation = as.numeric(par_sub$simulation)
    par_sub$truth = as.numeric(par_sub$truth)
    par_sub$est = as.numeric(par_sub$est)
    par_sub$var = as.numeric(par_sub$var)
    
    inds = which(ests$ve == sim_settings$ve[s] & 
                     ests$crossdur == sim_settings$crossdur[s])
    
    par_inds = which(pars$ve == sim_settings$ve[s] & 
                         pars$crossdur == sim_settings$crossdur[s])
    
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
}

# combine the results
res_ve = do.call(rbind, ve_list)
res_decay = do.call(rbind, decay_list)
res_pars = do.call(rbind, par_list)

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
    group_by(ve, crossdur, time, model) %>% 
    summarize(emp_var = var(est),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

decay_sum = 
    res_decay %>% 
    group_by(ve, crossdur, time, model) %>% 
    summarize(emp_var = var(est),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

par_sum = 
    res_pars %>% 
    group_by(ve, crossdur, model, param) %>% 
    summarize(emp_var = var(est), # sum((est - truth)^2) / (n()-1),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

# subset to get the true data generating mechanism and the pspline
ve_sum = 
    ve_sum %>% 
    filter(model == "pspline" | model == "linear") %>% 
    arrange(ve, crossdur, model, time) %>% 
    ungroup()

decay_sum = 
    decay_sum %>% 
    filter(model == "pspline" | model == "linear") %>% 
    arrange(ve, crossdur, model, time) %>% 
    ungroup()

ve_sum = 
    ve_sum[,c("ve", "crossdur", "model", "time", "emp_var", "covg")]

decay_sum = 
    decay_sum[,c("ve", "crossdur", "model", "time", "emp_var", "covg")]

res_sum = full_join(ve_sum, decay_sum, c("ve", "crossdur", "model", "time"))

# recode some of the variables
res_sum = 
    res_sum %>% 
    mutate(model = 
               case_when(model == "pspline" ~ "P-spline",
                         TRUE ~ "log-linear"),
           crossdur = 
               case_when(crossdur == 2 ~ "Two week crossover",
                         crossdur == 8 ~ "Two month crossover"))

par_sum = 
    par_sum %>% 
    mutate(model = 
               case_when(model == "pspline" ~ "P-spline",
                         TRUE ~ "log-linear"),
           crossdur = 
               case_when(crossdur == 2 ~ "Two week crossover",
                         crossdur == 8 ~ "Two month crossover"),
           param = 
               case_when(param == "intcpt" ~ "Intercept",
                         TRUE ~ "Linear trend"))

res_sum[,5:8] = round(res_sum[,5:8], digits = 3)
par_sum[,5:6] = round(par_sum[,5:6], digits = 3)

# pivot par_sum wider
par_sum = 
    par_sum %>% 
    pivot_wider(names_from = c("param"), values_from = c("emp_var", "covg"))
par_sum = par_sum[,c(1:3,4,6,5,7)]

# generate tables
res_tab <- 
    res_sum %>% 
    arrange(desc(ve), desc(crossdur), model) %>% 
    select(2:8)

res_tab$crossdur[-seq(1,40,by=10)] = " "
res_tab$model[-seq(1,40,by=5)] = " "

res_tab %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 2, "$VE(t)$" = 2, "$\\Delta VE(t)$" = 2)) %>%
    pack_rows("Vaccine efficacy constant at 75%", 1, 20) %>% 
    pack_rows("Vaccine efficacy wanes from 85% to 35% over 1.5 years", 21, 40)

# par table
partab <- 
    par_sum %>% 
    arrange(desc(ve), desc(crossdur)) %>% 
    pivot_wider(names_from = c("param"), values_from = c("emp_var", "covg"))


partab$crossdur[seq(2,8,by=2)] = " "

partab[,c(3,4,6,5,7)] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 1, 
                       "Intercept" = 2, "Linear Trend" = 2)) %>%
    pack_rows("Vaccine efficacy constant at 75%", 1, 4) %>% 
    pack_rows("Vaccine efficacy wanes from 85% to 35% over 1.5 years", 5,8) %>% 
    pack_rows("Two week crossover", 1, 2) %>% 
    pack_rows("Two month crossover", 3, 4) %>% 
    pack_rows("Two week crossover", 5, 6) %>% 
    pack_rows("Two month crossover", 7, 8)

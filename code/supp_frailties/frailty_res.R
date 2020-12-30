library(knitr)
library(kableExtra)
library(here)
library(tidyverse)
library(rsimsum)
setwd(here("jon_stuff/frailty_sim"))

sim_settings <- 
    expand.grid(cross = c("time", "parallel"),
                ve = c("waning", "constant"),
                event_rate = c(50, 300),
                frailty_sd = c(1, 2))

pars <- expand.grid(cross = c("time", "parallel"),
                    ve = c("waning", "constant"),
                    event_rate = c(50, 300),
                    frailty_sd = c(1, 2),
                    param = c("intcpt", "slope"),
                    model = c("const", "linear", "pspline"),
                    simulation = 1:1e4)

frailties <-
    expand.grid(cross = c("time", "parallel"),
                ve = c("waning", "constant"),
                event_rate = c(50, 300),
                frailty_sd = c(1, 2),
                simulation = 1:1e4,
                arm = c("plac", "trt"),
                stat = c("mean", "sd", "X10.", "X25.", "X50.", "X75.", "X90."))

pars_cross <- 
    expand.grid(cross = c("time", "parallel"),
                ve = c("waning", "constant"),
                event_rate = c(50, 300),
                frailty_sd = c(1, 2),
                model = c("const", "linear", "pspline"),
                simulation = 1:1e4)

ests <- expand.grid(cross = c("time", "parallel"),
                    ve = c("waning", "constant"),
                    event_rate = c(50, 300),
                    frailty_sd = c(1, 2),
                    time = c(0.5, 1, 1.5, 2),
                    model = c("const", "linear", "pspline"),
                    simulation = 1:1e4)

par_list <- vector("list", nrow(sim_settings))
ve_list <- vector("list", nrow(sim_settings))
decay_list <- vector("list", nrow(sim_settings))
t_cross_list <- vector("list", nrow(sim_settings))
n_cross_list <- vector("list", nrow(sim_settings)) 
frailty_list <- vector("list", nrow(sim_settings))

for(s in seq_len(nrow(sim_settings))) {
    
    filename =
        paste0("frailty_sim_", 
               sim_settings$ve[s], "_",
               sim_settings$cross[s], "_",
               sim_settings$event_rate[s], "_",
               sim_settings$frailty_sd[s], ".Rds")
    
    sim_res = readRDS(filename)
    
    # get n_cross and t_cross
    n_cross_list[[s]] <- sapply(sim_res, "[[", "n_cross")
    t_cross_list[[s]] <- sapply(sim_res, "[[", "t_cross")
    
    # frailty stats
    frailty_sub = 
        data.frame(expand.grid(
            stat = c("mean", "sd", "X10.", "X25.", "X50.", "X75.", "X90."),
            arm = c("plac", "trt"),
            simulation = 1:1e4),
            est = do.call(c, lapply(sim_res, function(x) unlist(x$frailty_sum))))
    
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
    
    inds = which(ests$cross == sim_settings$cross[s] & 
                     ests$ve == sim_settings$ve[s] & 
                     ests$event_rate == sim_settings$event_rate[s] & 
                     ests$frailty_sd == sim_settings$frailty_sd[s])
    
    par_inds = which(pars$cross == sim_settings$cross[s] & 
                         pars$ve == sim_settings$ve[s] & 
                         pars$event_rate == sim_settings$event_rate[s] & 
                         pars$frailty_sd == sim_settings$frailty_sd[s])
    
    frailty_inds = 
        which(frailties$cross == sim_settings$cross[s] & 
                  frailties$ve == sim_settings$ve[s] & 
                  frailties$event_rate == sim_settings$event_rate[s] & 
                  frailties$frailty_sd == sim_settings$frailty_sd[s])
    
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
    
    frailty_list[[s]] <- 
        left_join(frailties[frailty_inds,], frailty_sub,
                  by = c("stat", "arm", "simulation"))
}

# combine the results
res_ve = do.call(rbind, ve_list)
res_decay = do.call(rbind, decay_list)
res_pars = do.call(rbind, par_list)
res_frailies = do.call(rbind, frailty_list)

# summarize crossover times and event counts
cross_stats = 
    data.frame(
        sim_settings,
        n_cross_mean = sapply(n_cross_list, mean, na.rm = T),
        n_cross_sd   = sapply(n_cross_list, sd, na.rm = T),
        t_cross_mean = sapply(t_cross_list, mean, na.rm = T),
        t_cross_mean = sapply(t_cross_list, sd, na.rm = T)
    )

# frailty stats
frailty_stats = 
    res_frailies %>% 
    group_by(cross, ve, event_rate, frailty_sd, arm, stat) %>% 
    # subset(stat == "mean") %>% 
    summarize(mean = round(exp(mean(log(est))), digits = 3))

frailty_res = 
    frailty_stats %>% 
    pivot_wider(names_from = c("stat"),
                values_from = c("mean")) %>% 
    arrange(desc(ve), event_rate, frailty_sd, cross, arm) %>% 
    mutate(cross = as.character(cross),
           ve = as.character(ve),
           arm = as.character(arm)) %>% 
    mutate(cross = case_when(cross == "time" ~ "Cross at 1 year",
                             cross == "parallel" ~ "Parallel trial"),
           arm = case_when(arm == "plac" ~ "Original placebo arm",
                           arm == "trt" ~ "Original vaccine arm"),
           event_rate = case_when(event_rate == 50 ~ "Low baseline hazard",
                                  event_rate == 300 ~ "High baseline hazard"),
           frailty_sd = case_when(frailty_sd == 1 ~ "Low frailty variance",
                                  frailty_sd == 2 ~ "High frailty variance"))

frailty_res$cross[seq(2,nrow(frailty_res), by = 2)] = " "
frailty_res$event_rate[-seq(1,32,by=8)] = " "
frailty_res$frailty_sd[-seq(1,32,by=4)] = " "

frailty_res[,c(3,4,1,5,6,7,9:11)] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    pack_rows("VE constant at 75%", 1, 16) %>% 
    pack_rows("VE wanes from 85% to 35% over 1.5 years", 17, 32) 

# summarize
ve_sum = 
    res_ve %>% 
    group_by(cross, ve, event_rate, frailty_sd, time, model) %>% 
    summarize(bias = mean(est - truth),
              emp_var = var(est),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

decay_sum = 
    res_decay %>% 
    group_by(cross, ve, event_rate, frailty_sd, time, model) %>% 
    summarize(bias = mean(est - truth),
              emp_var = var(est),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

par_sum = 
    res_pars %>% 
    group_by(cross, ve, event_rate, frailty_sd, model, param) %>% 
    summarize(bias = mean(est - truth),
              emp_var = var(est), 
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var)))) 

par_sum = par_sum[complete.cases(par_sum),]

# subset to get the true data generating mechanism and the pspline
ve_sum = 
    ve_sum %>% 
    filter(model == "pspline" | model == "linear") %>% 
    arrange(ve, event_rate, frailty_sd, cross, model, time) %>% 
    ungroup()

decay_sum = 
    decay_sum %>% 
    filter(model == "pspline" | model == "linear") %>% 
    arrange(ve, event_rate, frailty_sd, cross, model, time) %>% 
    ungroup()

ve_sum = 
    ve_sum[,c("ve", "event_rate", "frailty_sd", 
              "cross", "model", "time", "bias")]

decay_sum = 
    decay_sum[,c("ve", "event_rate", "frailty_sd",
                 "cross", "model", "time", "bias")]

res_sum = full_join(ve_sum, decay_sum, c("ve", "event_rate", "frailty_sd", 
                                         "cross", "model", "time"))

# recode some of the variables
res_sum = 
    res_sum %>% 
    mutate(cross = 
               case_when(cross == "events" ~ "Cross at 150 cases",
                         cross == "time" ~ "Cross at 1 year",
                         TRUE ~ "Parallel trial"),
           model = 
               case_when(model == "pspline" ~ "P-spline",
                         model == "linear" ~ "log-linear",
                         model == "const" ~ "constant"))

par_sum$cross = as.character(par_sum$cross)
par_sum = 
    par_sum %>% 
    mutate(cross = 
               case_when(cross == "events" ~ "Cross at 150 cases",
                         cross == "time" ~ "Cross at 1 year",
                         cross == "parallel" ~ "Parallel trial",
                         TRUE ~ cross),
           model = 
               case_when(model == "pspline" ~ "P-spline",
                         model == "linear" ~ "log-linear",
                         model == "const" ~ "Constant VE"),
           param = 
               case_when(param == "intcpt" ~ "Intercept",
                         TRUE ~ "Linear trend"))

res_sum[,7:8] = round(res_sum[,7:8], digits = 3)
par_sum[,7:9] = round(par_sum[,7:9], digits = 3)

# pivot res_sum wider
res_sum = 
    res_sum %>% 
    pivot_wider(names_from = c("cross"), values_from = c("bias.x", "bias.y"))

res_sum = res_sum[,c(1:5, 6,8,7,9)]
names(res_sum)[6:9] = c("Cross_VE", "Cross_decay", "Parallel_VE", "Parallel_decay")

# pivot par_sum wider
par_sum = 
    par_sum[,1:7] %>% 
    pivot_wider(names_from = c("param"), values_from = c("bias"))

# generate tables
tab_const <- 
    res_sum %>% 
    filter(ve == "constant") %>%
    mutate(frailty_sd = case_when(frailty_sd == 1 ~ "Low",
                                  frailty_sd == 2 ~ "High"),
           event_rate = case_when(event_rate == 50 ~ "Low",
                                  event_rate == 300 ~ "High")) %>% 
    select(2:9)

tab_waning <- 
    res_sum %>% 
    filter(ve == "waning") %>%
    mutate(frailty_sd = case_when(frailty_sd == 1 ~ "Low",
                                  frailty_sd == 2 ~ "High"),
           event_rate = case_when(event_rate == 50 ~ "Low",
                                  event_rate == 300 ~ "High")) %>% 
    select(2:9)

tab_const$frailty_sd[-seq(1,nrow(tab_const), by = 8)] = " "
tab_const$model[-seq(1,nrow(tab_const), by = 4)] = " "

tab_waning$frailty_sd[-seq(1,nrow(tab_waning), by = 8)] = " "
tab_waning$model[-seq(1,nrow(tab_waning), by = 4)] = " "


tab_const[,-1] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 3, "Crossover" = 2, "Parallel" = 2)) %>%
    pack_rows("Low baseline hazard", 1, 16) %>% 
    pack_rows("High baseline hazard", 17, 32)


tab_waning[,-1] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 3, "Crossover" = 2, "Parallel" = 2)) %>%
    pack_rows("Low baseline hazard", 1, 16) %>% 
    pack_rows("High baseline hazard", 17, 32)



# par table
par_sum$cross = factor(par_sum$cross, 
                        levels = c("Cross at 1 year", "Parallel trial"))

partab <- 
    par_sum %>% 
    arrange(desc(ve), event_rate, frailty_sd, cross, model) %>% 
    mutate(frailty_sd = case_when(frailty_sd == 1 ~ "Low",
                                  frailty_sd == 2 ~ "High"),
           event_rate = case_when(event_rate == 50 ~ "Low",
                                  event_rate == 300 ~ "High"))

partab$cross = as.character(partab$cross)
partab$model = as.character(partab$model)

partab$cross[-seq(1,nrow(partab),by=3)] = " "
partab$event_rate[-seq(1,48,by=12)] = " "
partab$frailty_sd[-seq(1,48,by=6)] = " "

partab[,c(3,4,1,5:7)] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    pack_rows("Vaccine efficacy constant at 75%", 1, 24) %>% 
    pack_rows("Vaccine efficacy wanes from 85% to 35%", 25, 48)

library(knitr)
library(kableExtra)
library(here)
library(tidyverse)
library(rsimsum)
setwd(here("jon_stuff/pure_trials"))

pars <- expand.grid(design = c("uniform", "parallel", "cross1yr"),
                    param = c("intcpt", "slope"),
                    model = c("linear", "pspline"),
                    simulation = 1:1e4)

ests <- expand.grid(design = c("uniform", "parallel", "cross1yr"),
                    param = c("intcpt", "slope"),
                    model = c("linear", "pspline"),
                    time = c(0, 0.5, 1, 1.5, 2),
                    simulation = 1:1e4)

par_list <- vector("list", 3)
ve_list <- vector("list", 3)
decay_list <- vector("list", 3) 

for(s in seq_len(3)) {
    
    filename = paste0("pure_sim_", s, ".Rds")
    
    sim_res = readRDS(filename)
    
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
    
    inds = which(ests$design == 
                     case_when(s == 1 ~ "uniform",
                               s == 2 ~ "parallel",
                               s == 3 ~ "cross1yr"))
    par_inds = which(pars$design == 
                         case_when(s == 1 ~ "uniform",
                                   s == 2 ~ "parallel",
                                   s == 3 ~ "cross1yr"))
    
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

# summarize
ve_sum = 
    res_ve %>% 
    group_by(design, time, model) %>% 
    summarize(emp_var = var(est),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

decay_sum = 
    res_decay %>% 
    group_by(design, time, model) %>% 
    summarize(emp_var = var(est),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

par_sum = 
    res_pars %>% 
    group_by(design, model, param) %>% 
    summarize(emp_var = var(est), # sum((est - truth)^2) / (n()-1),
              covg = mean(truth > (est - 1.96 * sqrt(var)) & 
                              truth < (est + 1.96 * sqrt(var))))

# subset to get the true data generating mechanism and the pspline
ve_sum = 
    ve_sum %>% 
    arrange(design, model, time) %>% 
    ungroup()

decay_sum = 
    decay_sum %>% 
    arrange(design, model, time) %>% 
    ungroup()

ve_sum = 
    ve_sum[,c("design", "model", "time", "emp_var", "covg")]

decay_sum = 
    decay_sum[,c("design", "model", "time", "emp_var", "covg")]

res_sum = full_join(ve_sum, decay_sum, c("design", "model", "time"))

# recode some of the variables
res_sum = 
    res_sum %>% 
    mutate(design = 
               case_when(design == "uniform" ~ "Uniform crossover",
                         design == "cross1yr" ~ "Crossover at 1 year",
                         TRUE ~ "Parallel trial"),
           model = 
               case_when(model == "pspline" ~ "P-spline",
                         TRUE ~ "log-linear"))

par_sum = 
    par_sum %>% 
    mutate(design = 
               case_when(design == "uniform" ~ "Uniform crossover",
                         design == "cross1yr" ~ "Crossover at 1 year",
                         TRUE ~ "Parallel trial"),
           model = 
               case_when(model == "pspline" ~ "P-spline",
                         TRUE ~ "log-linear"),
           param = 
               case_when(param == "intcpt" ~ "Intercept",
                         TRUE ~ "Linear trend"))



res_sum[,4:7] = round(res_sum[,4:7], digits = 3)
par_sum[,4:5] = round(par_sum[,4:5], digits = 3)

# pivot par_sum wider
par_sum$design = factor(par_sum$design,
                        levels = c("Uniform crossover",
                                   "Crossover at 1 year", 
                                   "Parallel trial"))
par_sum = 
    par_sum %>% 
    pivot_wider(names_from = c("param"), values_from = c("emp_var", "covg"))
par_sum = par_sum[,c(1:2,3,5,4,6)]

# generate tables
res_sum$design[-seq(1,30,by=10)] = " "
res_sum$model[-seq(1,30,by=5)] = " "

res_sum[c(1:10,21:30,11:20),-1] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 2, "$VE(t)$" = 2, "$\\Delta VE(t)$" = 2)) %>%
    pack_rows("Continuous uniform crossover", 1, 10) %>% 
    pack_rows("Crossover at one year", 11, 20) %>% 
    pack_rows("Parallel trial", 21, 30) 

par_sum$design[seq(2,6,by=2)] = " "

par_sum[c(1:2,5:6,3:4),-1] %>% 
    knitr::kable(format = "latex", booktabs = T) %>% 
    kable_styling() %>% 
    pack_rows("Continuous uniform crossover", 1, 2) %>% 
    pack_rows("Crossover at one year", 3, 4) %>% 
    pack_rows("Parallel trial", 5, 6) 
    add_header_above(c(" " = 2, 
                       "log-linear model" = 2, "P-spline model" = 2))

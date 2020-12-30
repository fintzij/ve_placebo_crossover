library(simsurv) # simulates data
library(survival)
library(tidyverse)
library(ggthemes)
library(knitr)
library(patchwork)
library(kableExtra)

rep_est = function(est, se, digs = 2) {
    res = vector("list", length(est))
    est_r = unname(round(est, digits = digs))
    ci_r = unname(round(cbind(est - 1.96 * se, est + 1.96 * se), digits = digs))
    
    for(s in seq_along(est)) {
        res[[s]] = 
            paste0(est_r[s], " (95% CI: ", ci_r[s,1], ", ", ci_r[s,2], ")")
    }
    
    return(unlist(res))
}

# set up parameters
t_cross = 1
n_cross = NULL

# parameters - 100/6mo, 150/6mo, 100/6mo, 50/6mo in placebo
# and 50/6mo, 75/6mo, 50/6mo, 25/6mo in year 2
# h = function(r) (pexp(0.5, r) * 1.5e4 - 50)^2
# optimize(h, interval = c(0, 1))
sim_pars_waning = 
    c(VE_max = log(1 - 0.85),
      VE_decay = 0.977558,
      lambda_1 = 0.0134, # 100 / 6 months
      lambda_2 = 0.02,   # 150 / 6 months
      lambda_3 = 0.0134, # 100 / 6 months
      lambda_4 = 0.0067, # 50  / 6 months
      lambda_5 = 0.0067, # 50 / 6 months
      lambda_6 = 0.01,   # 75 / 6 months
      lambda_7 = 0.0067, # 50 / 6 months
      lambda_8 = 0.0033, # 25 / 6 mo
      lambda_9 = 0.0067) # 25 / 6 mo

sim_pars_const = 
    c(VE_max = log(1 - 0.75),
      VE_decay = 0,
      lambda_1 = 0.0134, # 100 / 6 months
      lambda_2 = 0.02,   # 150 / 6 months
      lambda_3 = 0.0134, # 100 / 6 months
      lambda_4 = 0.0067, # 50  / 6 months
      lambda_5 = 0.0067, # 50 / 6 months
      lambda_6 = 0.01,   # 75 / 6 months
      lambda_7 = 0.0067, # 50 / 6 months
      lambda_8 = 0.0033, # 25 / 6 mo
      lambda_9 = 0.0067) # 50 / 6 mo

# cumulative hazard in a time-homogeneous interval
c_haz_int = 
    function(lambda_t, t_l, t_r, trt,  
             VE_max = NULL, VE_decay = NULL, t_trt = NULL) {
        
        # cumulative hazard over any interval for unvaccinated, always
        if(trt == 0) {
            c_haz = lambda_t * (t_r - t_l)
        } else {
            if(VE_decay == 0) {
                c_haz = lambda_t * exp(VE_max) * (t_r - t_l)
            } else {
                c_haz = 
                    lambda_t * exp(VE_max) / VE_decay * 
                    (exp(VE_decay * (t_r - t_trt)) - exp(VE_decay * (t_l - t_trt)))
            }
        }
        return(c_haz)
    }

# cumulative hazard function
cumul_haz = function(t, x, betas, ...) {
    
    # initialize the hazard vector
    c_haz_t = rep(0.0, length(t))
    
    # get time when hazard starts accumulating
    tstart = max(betas[["t0"]], x[["entry"]])
    
    # baseline hazard intervals
    haz_ints = seq(0, 2.25, by=0.25)
    times = sort(c(unique(c(betas[["t0"]], x[["entry"]], x[["crossover"]],
                            haz_ints, t))))
    endpoints_l = head(times, length(times)-1)
    endpoints_r = tail(times, length(times)-1)
    
    # which intervals are contributing
    contribs = endpoints_l >= tstart & endpoints_r <= min(t, 2.25)
    
    # treatment indicators at left endpoints
    trt_inds = ifelse((x[["trt"]] == 1 & endpoints_l >= x[["entry"]]) |
                          (x[["trt"]] == 0 & endpoints_l >= x[["crossover"]]), 1, 0)
    
    # now compute the cumulative hazards
    for(s in seq_along(t)) {
        
        # initialize vector of cumulative hazards
        c_haz_ints_s = rep(0.0, length(endpoints_l))
        
        # compute the hazards for each interval
        for(i in seq_along(endpoints_l)) {
            if(contribs[i]) {
                
                # grab the baseline hazard
                lambda_t = 
                    betas[[paste0("lambda_", 
                                  findInterval(endpoints_l[i], 
                                               haz_ints, all.inside = T))]]
                
                # set time of treatment administration
                t_trt = ifelse(x[["trt"]] == 1, x[["entry"]], x[["crossover"]])
                
                # cumulative hazard contribution
                c_haz_ints_s[i] = 
                    c_haz_int(lambda_t = lambda_t, 
                              t_l = endpoints_l[i],
                              t_r = endpoints_r[i], 
                              trt = trt_inds[i],
                              VE_max = betas[["VE_max"]],
                              VE_decay = betas[["VE_decay"]],
                              t_trt = t_trt)
            }
        }
        
        c_haz_t[s] = sum(c_haz_ints_s)
    }
    
    return(c_haz_t)
}

sim_one <- function(n,
                    pars,
                    cumul_haz, 
                    n_cross,
                    t_cross,
                    enrolldur = 12/52, # enrollment period
                    crossdur = 4/52,  # crossover period
                    pvacc = 0.5) {
    
    if(!is.null(t_cross)) {
        
        # create a data frame with subject IDs and treatment indicators
        covars <- data.frame(id = seq_len(n),
                             trt = rbinom(n, 1, pvacc),
                             entry = runif(n, 0, enrolldur),
                             crossover = 0)
        
        # simulate crossover times for placebo group
        covars$crossover[covars$trt == 0] = 
            t_cross + runif(sum(covars$trt == 0), 1e-5, 1e-5 + crossdur)
        
        # simulate event times
        dat <- simsurv::simsurv(betas = c(pars, t0 = 0), 
                                x = covars,
                                cumhazard = cumul_haz,
                                maxt = 2.25,
                                rootfun = qlogis,
                                ids = covars$id,
                                idvar = "id")
        
        dat <- merge(covars, dat) %>% arrange(eventtime)
        
    } else if(!is.null(n_cross)) {
        
        # create a data frame with subject IDs and treatment indicators
        covars_1 <- data.frame(id = seq_len(n),
                               trt = rbinom(n, 1, pvacc),
                               entry = runif(n, 0, enrolldur),
                               crossover = Inf)
        
        # simulate event times
        dat_1 <- 
            simsurv::simsurv(betas = c(pars, t0 = 0),
                             x = covars_1,
                             cumhazard = cumul_haz,
                             maxt = 2.25,
                             rootfun = qlogis,
                             ids = covars_1$id,
                             idvar = "id")
        
        # sort by event times
        dat_1 <- dat_1 %>% arrange(eventtime)
        n_events = cumsum(dat_1$status)
        
        # grab the crossover time
        t_cross = dat_1$eventtime[n_events == n_cross][1]
        
        # cache the pre-crossover data
        dat_1 <- dat_1 %>% filter(dat_1$eventtime <= t_cross)
        dat_sim_1 = merge(covars_1, dat_1) %>%  arrange(eventtime)
        
        # subset to get covariates for next part of the simulation
        covars_2 = covars_1[!covars_1$id %in% dat_sim_1$id, ]
        
        # generate crossover times
        covars_2$crossover[covars_2$trt == 0] = 
            t_cross + runif(sum(covars_2$trt == 0), 1e-5, 1e-5 + crossdur)
        
        # keep simulating
        dat_2 <- 
            simsurv::simsurv(betas = c(pars, t0 = t_cross),
                             x = covars_2,
                             cumhazard = cumul_haz,
                             maxt = 2.25,
                             rootfun = qlogis,
                             ids = covars_2$id,
                             idvar = "id")
        
        # sort by event time
        dat_2 <- dat_2 %>% arrange(id)
        
        # merge covars_2 and dat_2
        dat_sim_2 = merge(covars_2, dat_2) %>% arrange(eventtime)
        
        # combine simulated datasets
        dat = rbind(dat_sim_1, dat_sim_2)
        dat = dat %>% arrange(eventtime)  
    }
    
    dat$id = seq_along(dat$id)
    
    # censor people at two years post-study entry
    cens_inds = dat$eventtime > (dat$entry + 2)
    dat$eventtime[cens_inds] = dat$entry[cens_inds] + 2
    dat$status[cens_inds] = 0
    
    return(dat)
}

# to censor people at crossover
cens_cross <- function(dat, n_cross, t_cross) {
    
    if(!is.null(n_cross)) {
        
        # get the crossover time
        t_cross = dat[n_cross, "eventtime"]
        
        # censor events after crossover
        dat$status[(n_cross + 1):nrow(dat)] = 0
        dat$eventtime[(n_cross + 1):nrow(dat)] = t_cross
        
        # get rid of the crossover indicator
        dat$crossover = NULL  
        
    } else {
        
        # get the crossover index
        cross_ind = which(dat[,"eventtime"] > t_cross)[1]
        
        # censor events after crossover
        dat$status[cross_ind:nrow(dat)] = 0
        dat$eventtime[cross_ind:nrow(dat)] = t_cross
        
        # get rid of the crossover indicator
        dat$crossover = NULL
    }
    
    # drop records where people were censored prior to entry
    dat <- subset(dat, eventtime > entry)
    
    # censor people at two years post-study entry
    cens_inds = dat$eventtime > (dat$entry + 2)
    dat$eventtime[cens_inds] = dat$entry[cens_inds] + 2
    dat$status[cens_inds] = 0
    
    # rename columns = 
    names(dat)[1:5] = c("id", "vacc", "tstart", "tstop", "status")
    
    return(dat)
}

reshape_dat <- function(dat, align_entry = FALSE) {
    
    # align times so that entry times correspond to time zero for each participant
    if(align_entry) {
        dat <- 
            dat %>% 
            group_by(id) %>% 
            mutate(crossover = ifelse(trt == 1, 0, crossover - entry),
                   eventtime = eventtime - entry) %>% 
            mutate(entry = 0) %>% 
            ungroup() %>% 
            as.data.frame()
    }
    
    # two rows for each placebo, one for each in the treatment arm
    dat_new <- 
        data.frame(id = unlist(mapply(rep, x = dat$id, 
                                      each = ifelse(dat$trt == 1, 1, 2))),
                   tstart = unlist(lapply(seq_len(nrow(dat)), 
                                          function(x) {
                                              if(dat$trt[x] == 1) {
                                                  dat$entry[x]  
                                              } else {
                                                  c(dat$entry[x], dat$crossover[x])
                                              }})),
                   tstop = NA,
                   status = NA,
                   assg = 0,
                   vacc = 0,
                   tvacc = 0)
    
    # fill out the interval endpoints and status
    tmax = max(dat$eventtime)
    for(s in seq_along(unique(dat_new$id))) {
        
        inds <- which(dat_new$id == dat$id[s])
        
        if(dat$trt[s] == 1) {
            
            dat_new$assg[inds] = 1
            dat_new$vacc[inds] = 1
            dat_new$tvacc[inds] = dat$entry[s]
            dat_new$tstop[inds] = dat$eventtime[s]
            dat_new$status[inds] = dat$status[s]
            
        } else {
            
            dat_new$assg[inds] = c(0,0)
            dat_new$vacc[inds] = c(0,1)
            dat_new$tvacc[inds] = c(dat$crossover[s], dat$crossover[s])
            
            if(dat$eventtime[s] < dat$crossover[s]) {
                
                dat_new$tstop[inds]  = c(dat$eventtime[s], NA)
                dat_new$status[inds] = c(dat$status[s], NA)
                
            } else {
                
                dat_new$tstop[inds] = c(dat$crossover[s], dat$eventtime[s])
                dat_new$status[inds] = c(0, dat$status[s])
            }
        }
    }
    
    dat_new = 
        dat_new[complete.cases(dat_new),] %>% 
        arrange(id, tstop)
    
    return(dat_new)
}

# Simulate trials -----------------------------------------------------------------------------

# set.seed
set.seed(1834)

# simulate trial
sim_waning <- 
    sim_one(n = 3e4,
            pars = sim_pars_waning, 
            n_cross = n_cross,
            t_cross = t_cross,
            cumul_haz = cumul_haz)

sim_const <- 
    sim_one(n = 3e4,
            pars = sim_pars_const, 
            n_cross = n_cross,
            t_cross = t_cross,
            cumul_haz = cumul_haz)

# reshape data
simr_const <- reshape_dat(sim_const, align_entry = FALSE)
simr_waning <- reshape_dat(sim_waning, align_entry = FALSE)

simr_const_aligned <- reshape_dat(sim_const, align_entry = TRUE)
simr_waning_aligned <- reshape_dat(sim_waning, align_entry = TRUE)

# simulation table
sim_settings = 
    expand.grid(ve = c("constant", "waning"),
                timepoint = c("interim", "cross", "end"),
                model = c("lin", "ps"))

par_res = 
    data.frame(
        expand.grid(param = c("intercept", "slope", "lrt"),
                    model = c("lin", "ps"),
                    timepoint = c("interim", "cross", "end"),
                    ve = c("constant", "waning")),
        est = 0)

simstats = 
    data.frame(t_look = ceiling(c(sim_const$eventtime[150],
                                  sim_waning$eventtime[150])*365.25),
               n_look = 
                   c(with(sim_const[1:150,], 
                          paste0(sum(status[trt == 0]), "_", sum(status[trt == 1]))),
                     with(sim_waning[1:150,], 
                          paste0(sum(status[trt == 0]), "_", sum(status[trt == 1])))),
               n_cross = 
                   c(with(subset(sim_const, eventtime <= 1), 
                          paste0(sum(status[trt == 0]), "_", sum(status[trt == 1]))),
                     with(subset(sim_waning, eventtime <= 1), 
                          paste0(sum(status[trt == 0]), "_", sum(status[trt == 1])))),
               n_tot = 
                   c(with(sim_const, 
                          paste0(sum(status[trt == 0]), "_", sum(status[trt == 1]))),
                     with(sim_waning, 
                          paste0(sum(status[trt == 0]), "_", sum(status[trt == 1])))), 
               row.names = c("constant", "waning"))

# fill out the results
for(s in seq_len(nrow(sim_settings))) {
    
    if(sim_settings$timepoint[s] == "end") {
            if(sim_settings$ve[s] == "constant") {
                dat = simr_const
            } else {
                dat = simr_waning
            }
    } else if(sim_settings$timepoint[s] == "interim") {
            if(sim_settings$ve[s] == "constant") {
                dat = cens_cross(sim_const, n_cross = 150, t_cross = NULL)
            } else {
                dat = cens_cross(sim_waning, n_cross = 150, t_cross = NULL)
            }
    } else if(sim_settings$timepoint[s] == "cross") {
            if(sim_settings$ve[s] == "constant") {
                dat = cens_cross(sim_const, n_cross = NULL, t_cross = 1)
            } else {
                dat = cens_cross(sim_waning, n_cross = NULL, t_cross = 1)
            }
    }
    
    # add tvacc to censored data
    if(sim_settings$timepoint[s] != "end") 
        dat$tvacc = ifelse(dat$vacc == 0, Inf, dat$tstart)    
    
    # fit model
    mod = if(sim_settings$model[s] == "ps") {
        coxph(formula =
                  Surv(time = tstart, 
                       time2 = tstop, 
                       event = status) ~ vacc + tt(tvacc),
              tt = function(tvacc, t, ...) {
                  pspline(pmax(0, t - tvacc), df = 3, nterm = 8)
              },
              data = dat)
    } else {
        coxph(formula =
                  Surv(time = tstart, 
                       time2 = tstop, 
                       event = status) ~ vacc + tt(tvacc),
              tt = function(tvacc, t, ...) {
                  pmax(0, t - tvacc)
              },
              data = dat)
    }
    
    # model without waning VE
    mod_nowaning = update(mod, .~. - tt(tvacc))
    
    # likelihood ratio test for waning VE
    lrt_pval = 
        round(1 - pchisq(2 * (mod$loglik - mod_nowaning$loglik)[2], 
                   df = ifelse(sim_settings$model[s] == "lin", 1,
                               mod$df[2] - 1)), digits = 5)
    
    # estimates and CIs
    ests = rep_est(summary(mod)$coefficients[1:2,"coef"], 
                   summary(mod)$coefficients[1:2,"se(coef)"])
    
    # save estimates and LRT p-value
    par_res$est[which(par_res$ve == sim_settings$ve[s] & 
                          par_res$timepoint == sim_settings$timepoint[s] & 
                          par_res$model == sim_settings$model[s])] = c(ests, lrt_pval)
    
}

res = cbind(par_res[1:18,-c(4)], par_res[19:36, -c(1:4)])[,c(3:1,4,5)]
res = res %>% 
    mutate(timepoint = 
               case_when(timepoint == "interim" ~ "Estimates at interim look", 
                         timepoint == "cross" ~ "Estimates at 1 year crossover", 
                         timepoint == "end" ~ "Estimates at 2 year followup"),
           model = case_when(model == "lin" ~ "log-linear model", 
                             model == "ps" ~ "P-spline model"),
           param = 
               case_when(param == "intercept" ~ "Intercept",
                         param == "slope" ~ "Linear trend",
                         param == "lrt" ~ "LRT for time-varying VE"))
    
res$timepoint[-seq(1,36,by=6)] = " "
res$model[-seq(1,36,by=3)] = " "

res[,-c(1:2)] %>%
    kable(format = "latex", booktabs = T) %>% 
    pack_rows("Estimates at interim look", 1, 6) %>% 
    pack_rows("Estimates at 1 year crossover", 7, 12) %>% 
    pack_rows("Estimates at 2 year followup", 13, 18) %>% 
    pack_rows("log-linear model", 1, 3) %>% 
    pack_rows("P-spline model", 4, 6) %>% 
    pack_rows("log-linear model", 7, 9) %>% 
    pack_rows("P-spline model", 10, 12) %>% 
    pack_rows("log-linear model", 13, 15) %>% 
    pack_rows("P-spline model", 16, 18)


# Severe cases --------------------------------------------------------------------------------

set.seed(1834)
sim_const$severe = ifelse(sim_const$status == 0, 0, rbinom(nrow(sim_const), 1, 1/8))
sim_waning$severe = ifelse(sim_waning$status == 0, 0, rbinom(nrow(sim_waning), 1, 1/8))
simr_const$severe = 0
simr_waning$severe = 0

for(s in seq_along(sim_const$id)) {
    if(sim_const$severe[s] == 1) {
        inds = which(simr_const$id == sim_const$id[s])
        simr_const$severe[tail(inds, 1)] = 1
    }
    
    if(sim_waning$severe[s] == 1) {
        inds = which(simr_waning$id == sim_waning$id[s])
        simr_waning$severe[tail(inds, 1)] = 1
    }
}

par_res = 
    data.frame(
        expand.grid(param = c("intercept", "slope", "lrt"),
                    model = c("lin", "ps"),
                    timepoint = c("interim", "cross", "end"),
                    ve = c("constant", "waning")),
        est = 0)

simstats = 
    data.frame(t_look = ceiling(c(sim_const$eventtime[150],
                                  sim_waning$eventtime[150])*365.25),
               n_look = 
                   c(with(sim_const[1:150,], 
                          paste0(sum(severe[trt == 0]), "_", sum(severe[trt == 1]))),
                     with(sim_waning[1:150,], 
                          paste0(sum(severe[trt == 0]), "_", sum(severe[trt == 1])))),
               n_cross = 
                   c(with(subset(sim_const, eventtime <= 1), 
                          paste0(sum(severe[trt == 0]), "_", sum(severe[trt == 1]))),
                     with(subset(sim_waning, eventtime <= 1), 
                          paste0(sum(severe[trt == 0]), "_", sum(severe[trt == 1])))),
               n_tot = 
                   c(with(sim_const, 
                          paste0(sum(severe[trt == 0]), "_", sum(severe[trt == 1]))),
                     with(sim_waning, 
                          paste0(sum(severe[trt == 0]), "_", sum(severe[trt == 1])))), 
               row.names = c("constant", "waning"))

# fill out the results
for(s in seq_len(nrow(sim_settings))) {
    
    if(sim_settings$timepoint[s] == "end") {
        if(sim_settings$ve[s] == "constant") {
            dat = simr_const
        } else {
            dat = simr_waning
        }
    } else if(sim_settings$timepoint[s] == "interim") {
        if(sim_settings$ve[s] == "constant") {
            dat = cens_cross(sim_const, n_cross = 150, t_cross = NULL)
        } else {
            dat = cens_cross(sim_waning, n_cross = 150, t_cross = NULL)
        }
    } else if(sim_settings$timepoint[s] == "cross") {
        if(sim_settings$ve[s] == "constant") {
            dat = cens_cross(sim_const, n_cross = NULL, t_cross = 1)
        } else {
            dat = cens_cross(sim_waning, n_cross = NULL, t_cross = 1)
        }
    }
    
    # add tvacc to censored data
    if(sim_settings$timepoint[s] != "end") 
        dat$tvacc = ifelse(dat$vacc == 0, Inf, dat$tstart)    
    
    # fit model
    mod = if(sim_settings$model[s] == "ps") {
        coxph(formula =
                  Surv(time = tstart, 
                       time2 = tstop, 
                       event = severe) ~ vacc + tt(tvacc),
              tt = function(tvacc, t, ...) {
                  pspline(pmax(0, t - tvacc), df = 3, nterm = 8)
              },
              data = dat)
    } else {
        coxph(formula =
                  Surv(time = tstart, 
                       time2 = tstop, 
                       event = severe) ~ vacc + tt(tvacc),
              tt = function(tvacc, t, ...) {
                  pmax(0, t - tvacc)
              },
              data = dat)
    }
    
    # model without waning VE
    mod_nowaning = update(mod, .~. - tt(tvacc))
    
    # likelihood ratio test for waning VE
    lrt_pval = 
        round(1 - pchisq(2 * (mod$loglik - mod_nowaning$loglik)[2], 
                         df = ifelse(sim_settings$model[s] == "lin", 1,
                                     mod$df[2] - 1)), digits = 5)
    
    # estimates and CIs
    ests = rep_est(summary(mod)$coefficients[1:2,"coef"], 
                   summary(mod)$coefficients[1:2,"se(coef)"])
    
    # save estimates and LRT p-value
    par_res$est[which(par_res$ve == sim_settings$ve[s] & 
                          par_res$timepoint == sim_settings$timepoint[s] & 
                          par_res$model == sim_settings$model[s])] = c(ests, lrt_pval)
    
}

res = cbind(par_res[1:18,-c(4)], par_res[19:36, -c(1:4)])[,c(3:1,4,5)]
res = res %>% 
    mutate(timepoint = 
               case_when(timepoint == "interim" ~ "Estimates at interim look", 
                         timepoint == "cross" ~ "Estimates at 1 year crossover", 
                         timepoint == "end" ~ "Estimates at 2 year followup"),
           model = case_when(model == "lin" ~ "log-linear model", 
                             model == "ps" ~ "P-spline model"),
           param = 
               case_when(param == "intercept" ~ "Intercept",
                         param == "slope" ~ "Linear trend",
                         param == "lrt" ~ "LRT for time-varying VE"))

res$timepoint[-seq(1,36,by=6)] = " "
res$model[-seq(1,36,by=3)] = " "

res[,-c(1:2)] %>%
    kable(format = "latex", booktabs = T) %>% 
    pack_rows("Estimates at interim look", 1, 6) %>% 
    pack_rows("Estimates at 1 year crossover", 7, 12) %>% 
    pack_rows("Estimates at 2 year followup", 13, 18) %>% 
    pack_rows("log-linear model", 1, 3) %>% 
    pack_rows("P-spline model", 4, 6) %>% 
    pack_rows("log-linear model", 7, 9) %>% 
    pack_rows("P-spline model", 10, 12) %>% 
    pack_rows("log-linear model", 13, 15) %>% 
    pack_rows("P-spline model", 16, 18)


# Plot VE(t) ----------------------------------------------------------------------------------

fit_ps_const = 
  coxph(formula =
        Surv(time = tstart, 
             time2 = tstop, 
             event = status) ~ vacc + tt(tvacc),
      tt = function(tvacc, t, ...) {
        pspline(pmax(0, t - tvacc), df = 3, nterm = 8)
      },
      data = simr_const)

fit_lin_const = 
  coxph(formula =
          Surv(time = tstart, 
               time2 = tstop, 
               event = status) ~ vacc + tt(tvacc),
        tt = function(tvacc, t, ...) {
          pmax(0, t - tvacc)
        },
        data = simr_const)

fit_ps_waning = 
  coxph(formula =
          Surv(time = tstart, 
               time2 = tstop, 
               event = status) ~ vacc + tt(tvacc),
        tt = function(tvacc, t, ...) {
          pspline(pmax(0, t - tvacc), df = 3, nterm = 8)
        },
        data = simr_waning)

fit_lin_waning = 
  coxph(formula =
        Surv(time = tstart, 
             time2 = tstop, 
             event = status) ~ vacc + tt(tvacc),
      tt = function(tvacc, t, ...) {
        pmax(0, t - tvacc)
      },
      data = simr_waning)

# estimates of VE(t)
VE_ests = 
  expand.grid(time = seq(0,2,by=0.01),
              Model = c("log-linear", "P-spline"),
              Scenario = c("Constant VE", "Waning VE"))

VE_ests$truth = 
  ifelse(VE_ests$Scenario == "Constant VE", log(0.25),
         log(0.15) + 0.977558 * VE_ests$time)

times_lin = cbind(1, seq(0,2,by=0.01))
times_ps = cbind(1, pspline(seq(0,2,by=0.01), df = 3, nterm = 8))
centvec = rep(0, ncol(times_ps)); centvec[2] = 1 # for centering the spline
centvec = rep(1, ncol(times_ps)) %*% diag(centvec)

VE_ests$est = 
  c(times_lin %*% coef(fit_lin_const),
    apply(times_ps, 1, function(x) (x - centvec) %*% coef(fit_ps_const)),
    times_lin %*% coef(fit_lin_waning),
    apply(times_ps, 1, function(x) (x - centvec) %*% coef(fit_ps_waning)))

VE_ests$lower = 
  VE_ests$est - 1.96 * c(
    sqrt(apply(times_lin, 1, function(x) t(x) %*% vcov(fit_lin_const) %*% x)),
    sqrt(apply(times_ps, 1, function(x) (x - centvec) %*% vcov(fit_ps_const) %*% t(x - centvec))),
    sqrt(apply(times_lin, 1, function(x) t(x) %*% vcov(fit_lin_waning) %*% x)),
    sqrt(apply(times_ps, 1, function(x) (x - centvec) %*% vcov(fit_ps_waning) %*% t(x - centvec))))

VE_ests$upper = 
  VE_ests$est + 1.96 * c(
    sqrt(apply(times_lin, 1, function(x) t(x) %*% vcov(fit_lin_const) %*% x)),
    sqrt(apply(times_ps, 1, function(x) (x - centvec) %*% vcov(fit_ps_const) %*% t(x - centvec))),
    sqrt(apply(times_lin, 1, function(x) t(x) %*% vcov(fit_lin_waning) %*% x)),
    sqrt(apply(times_ps, 1, function(x) (x - centvec) %*% vcov(fit_ps_waning) %*% t(x - centvec))))

require(ggthemes)
VE_ests %>% 
  ggplot(aes(x = time)) + 
  geom_line(aes(y = truth, linetype = "Truth"), linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model), alpha = 0.4) + 
  geom_line(aes(y = est, colour = Model)) + 
  facet_grid(Scenario ~ Model, scales = "free") + 
  theme_minimal() +
  # scale_color_colorblind() +
  # scale_fill_colorblind() +
  scale_color_tableau(palette = "Tableau 10", direction = -1) +
  scale_fill_tableau(palette = "Tableau 10", direction = -1) +
  scale_y_continuous(breaks = seq(-3, 2.5, by = 1),
                     labels = round(1 - exp(seq(-3,2.5,by=1)), digits = 1)) +
  labs(x = "Years since vaccination", y = "Vaccine efficacy", 
       title = "Vaccine efficacy vs. years since vaccination") + 
  theme(legend.position = "none",
        text = element_text(family = "serif", size = 20)) -> VE_plot

# table with events
simr_events_aligned = 
  rbind(cbind(ve = "Constant VE",
              simr_const_aligned[simr_const_aligned$status == 1,]),
        cbind(ve = "Waning VE",
              simr_waning_aligned[simr_waning_aligned$status == 1,]))

simr_events_aligned <-
  simr_events_aligned %>% 
  mutate(Arm = case_when(assg == 0 & vacc == 0 ~ "Placebo",
                         assg == 0 & vacc == 1 ~ "Delayed vaccination",
                         assg == 1 & vacc == 1 ~ "Immediate vaccination")) %>% 
  as.data.frame()

simr_events_aligned$Arm = 
  factor(simr_events_aligned$Arm, levels = c("Placebo", "Delayed vaccination",
                                             "Immediate vaccination"))

simr_events_lumped = 
  data.frame(expand.grid(ve = c("Constant VE", "Waning VE"),
                         Arm = ))
  simr_events_aligned %>% 
  group_by(ve, Arm) %>% 
  summarise(count = table(findInterval(tstop, seq(0,2,by=0.25))))

simr_events_aligned %>% 
  ggplot() + 
  geom_histogram(aes(x = tstop, group = Arm, fill = Arm), 
           stat = "bin", binwidth = 0.25, 
           position = position_dodge(width = 0.15),
           alpha = 0.5,
           boundary = 0, 
           closed = "right") + 
  scale_fill_colorblind() + 
  scale_x_continuous(breaks = seq(0,2,by=0.25)) + 
  facet_grid(ve ~ .) + 
  labs(x = "Years since vaccination", y = "Number of events",
       title = "Number of events by treatment arm") + 
  theme_minimal() + 
  theme(text = element_text(family = "serif", size = 20),
        legend.position = "bottom",
        panel.grid.minor.x = element_blank()) -> event_plot

example_plots = event_plot / VE_plot 

ggsave("example_plots.pdf", example_plots, height = 10, width = 10)

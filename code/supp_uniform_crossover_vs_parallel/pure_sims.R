library(simsurv) # simulates data
library(survival)
library(tidyverse)
library(foreach)
library(doRNG)
library(doMPI)
library(doParallel)
library(future)

args <- commandArgs(TRUE)
print(args)
simnum <- as.numeric(args[1])

# make cluster
cl <- doMPI::startMPIcluster()
registerDoMPI(cl)

sim_settings <- c("uniform", "parallel", "cross1yr")

t_cross = 
  if(simnum == 1) {
    0
  } else if(simnum == 2) {
    Inf 
  } else {
    1
  }

n_cross = NULL

# parameters - 100/6mo, 150/6mo, 100/6mo, 50/6mo in placebo
# or 50/6mo, 75/6mo, 50/6mo, 25/6mo in year 2
h = function(r) (pexp(0.5, r) * 1.5e3 - 50)^2
optimize(h, interval = c(0, 1))

sim_pars = c(VE_max = log(1 - 0.75), VE_decay = 0, 
             lambda_1 = 0.138,
             lambda_2 = 0.138,
             lambda_3 = 0.138,
             lambda_4 = 0.138,
             lambda_5 = 0.138, 
             lambda_6 = 0.138,
             lambda_7 = 0.138,
             lambda_8 = 0.138)

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
  haz_ints = seq(0, 2, by=0.25)
  times = sort(c(unique(c(betas[["t0"]], x[["entry"]], x[["crossover"]],
                          haz_ints, t))))
  endpoints_l = head(times, length(times)-1)
  endpoints_r = tail(times, length(times)-1)
  
  # which intervals are contributing
  contribs = endpoints_l >= tstart & endpoints_r <= min(t, 2)
  
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
                    enrolldur = 0, # enrollment period
                    crossdur = ifelse(simnum == 3, 0, 2),  # crossover period
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
                            maxt = 2,
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

# Test cumulative hazard function -------------------------------------------------------------

# hazard function
# haz = function(t, x, betas, ...) {
# 
#   # initialize the hazard vector
#   haz_t = rep(0.0, length(t))
# 
#   for(s in seq_along(t)) {
# 
#     if(t[s] <= betas[["t0"]] | t[s] <= x[["entry"]] | t[s] > 2) {
#       haz_t[s] = 0
# 
#     } else {
# 
#       # number of periods elapsed
#       pds = min(t[s] %/% 0.25, 7)
#       lambda_t = betas[paste0("lambda_", pds + 1)]
# 
#       if(x[["trt"]] == 1) {
# 
#         haz_t[s] = lambda_t * exp(betas[["VE_max"]] + betas[["VE_decay"]] * (t[s] - x[["entry"]]))
# 
#       } else {
#         if(t[s] <= x[["crossover"]]) {
#           haz_t[s] = lambda_t
#         } else {
#           haz_t[s] = lambda_t * exp(betas[["VE_max"]] +
#                                       betas[["VE_decay"]] * (t[s] - x[["crossover"]]))
#         }
#       }
#     }
#   }
# 
#   return(haz_t)
# }
# 
# haz_check =
#   data.frame(
#     expand.grid(trt = c(0,1),
#                 entry = 0.1,
#                 crossover = c(0.1, 0.6),
#                 times = c(seq(0.05, 2.15, by = 0.2), 2, 2 + 1e-9)),
#     c_haz_fn = 0,
#     c_haz_man = 0)
# 
# for(s in seq_len(nrow(haz_check))) {
# 
#   x = c(trt = haz_check$trt[s],
#         entry = haz_check$entry[s],
#         crossover = haz_check$crossover[s])
# 
#   haz_check$c_haz_man[s] =
#     integrate(haz,
#               lower = 0,
#               upper = haz_check$times[s],
#               x = x,
#               betas = c(sim_pars, t0 = 0),
#               rel.tol = .Machine$double.eps^(2/3),
#               subdivisions = 1e3)$value
# 
#   haz_check$c_haz_fn[s] =
#     cumul_haz(t = haz_check$times[s], x = x, betas = c(sim_pars, t0 = 0))
# }
# haz_check$diff <- haz_check$c_haz_fn - haz_check$c_haz_man
# all(abs(haz_check$diff) < 1e-8) # cool, congrats to the team

# Simulate trials -----------------------------------------------------------------------------

# set.seed
set.seed(12511)    
iters = 1e4

# simulate trials
trials = 
  foreach(i = seq_len(iters),
          .export = ls(),
          .packages = c("simsurv", "survival")) %dorng% {
            
            # simulate trial
            sim <- 
              sim_one(n = 3e3, 
                      pars = sim_pars, 
                      n_cross = n_cross,
                      t_cross = t_cross,
                      cumul_haz = cumul_haz)
            
            # reshape data
            simr <- reshape_dat(sim, align_entry = F)
            
            # possible to get a zero length time interval, rule it out
            while(min(simr$tstop - simr$tstart) < 1e-6) {
              sim <- 
                sim_one(n = 3e3, 
                        pars = sim_pars, 
                        n_cross = n_cross,
                        t_cross = t_cross,
                        cumul_haz = cumul_haz)
              
              # reshape data
              simr <- reshape_dat(sim, align_entry = F)  
            }
            
            # get the number at crossover or time of crossover
            if(is.null(n_cross)) n_cross = which(sim$eventtime > t_cross)[1] - 1
            if(is.null(t_cross)) t_cross = sim$eventtime[n_cross]
            
            # get the number of events at times 0.5, 1, 1.5, and 2
            event_counts <- 
              data.frame(time = c(0.5, 1, 1.5, 2),
                         events = 
                           sapply(c(0.5, 1, 1.5, 2), 
                                  function(x) which(sim$eventtime >= x)[1] - 1))
            
            ### fit models
            sim_cens = NULL
            fit_cross = NULL

            # constant VE, whole trial
            fit_const <- 
              coxph(formula = 
                      Surv(time = tstart, 
                           time2 = tstop, event = status) ~ vacc,
                    data = simr)
            
            # linear model
            fit_lin <- 
              coxph(formula = 
                      Surv(time = tstart, 
                           time2 = tstop, 
                           event = status) ~ vacc + tt(tvacc),
                    tt = function(tvacc, t, ...) {
                      pmax(0, t - tvacc)
                    },
                    data = simr)
            
            # quadratic model
            fit_quad <- 
              coxph(formula = 
                      Surv(time = tstart, 
                           time2 = tstop, 
                           event = status) ~ vacc + tt(tvacc),
                    tt = function(tvacc, t, ...) {
                      cbind(tvacc1 = pmax(0, t - tvacc), 
                            tvacc2 = pmax(0, t - tvacc)^2)
                    },
                    data = simr)
            
            # p-spline model
            fit_ps <-
              coxph(formula =
                      Surv(time = tstart, 
                           time2 = tstop, 
                           event = status) ~ vacc + tt(tvacc),
                    tt = function(tvacc, t, ...) {
                      pspline(pmax(0, t - tvacc), df = 3, nterm = 8)
                    },
                    data = simr)
            
            ### parameter estimates
            par_ests = 
              data.frame(
                simulation = i,
                truth = c(log(0.25), log(0.25), 0, log(0.25), 0),
                model = c("const", "linear", "linear", "pspline", "pspline"),
                est = c(coef(fit_const), coef(fit_lin), 
                        summary(fit_ps)$coefficients[1:2,"coef"]),
                var = c(vcov(fit_const), diag(vcov(fit_lin)),
                        summary(fit_ps)$coefficients[1:2,"se(coef)"]^2),
                row.names = NULL)
            
            ### VE estimates
            VE_ests = 
              data.frame(
                simulation = i,
                model = rep(c("const", "linear", "quadratic", "pspline"), 
                            each = 6),
                time = c(0, 0.25, 0.5, 1, 1.5, 2),
                truth = sim_pars["VE_max"],
                est = 0, 
                var = 0, 
                row.names = NULL)
            
            VE_ests$est[VE_ests$model == "const"] = coef(fit_const)
            VE_ests$var[VE_ests$model == "const"] = vcov(fit_const)[1,1]
            
            # linear model
            times_lin = cbind(1, c(0, 0.25, 0.5, 1, 1.5, 2))
            VE_ests$est[VE_ests$model == "linear"] = times_lin %*% coef(fit_lin)
            VE_ests$var[VE_ests$model == "linear"] = 
              apply(times_lin, 1, function(x) t(x) %*% vcov(fit_lin) %*% x)
            
            # quadratic model
            times_quad = 
              cbind(1, c(0, 0.25, 0.5, 1, 1.5, 2),
                    c(0, 0.25, 0.5, 1, 1.5, 2)^2)
            VE_ests$est[VE_ests$model == "quadratic"] = 
              times_quad %*% coef(fit_quad)
            VE_ests$var[VE_ests$model == "quadratic"] = 
              apply(times_quad, 1, function(x) t(x) %*% vcov(fit_quad) %*% x)
            
            # pspline model
            times_ps =
              cbind(1, pspline(c(0, 0.25, 0.5, 1, 1.5, 2), df = 3, nterm = 8))
            centvec = rep(0, ncol(times_ps)); centvec[2] = 1 # for centering the spline
            centvec = rep(1, ncol(times_ps)) %*% diag(centvec)
            
            VE_ests$est[VE_ests$model == "pspline"] = 
              apply(times_ps, 1, function(x) (x - centvec) %*% coef(fit_ps))
            VE_ests$var[VE_ests$model == "pspline"] = 
              apply(times_ps, 1, 
                    function(x) (x - centvec) %*% vcov(fit_ps) %*% t(x - centvec))
            
            ### VE decays
            VE_decays = 
              data.frame(
                simulation = i,
                model = rep(c("linear", "quadratic", "pspline"), each = 4),
                time = c(0.5, 1, 1.5, 2.0),
                truth = 0,
                est = 0,
                var = 0,
                row.names = NULL)
            
            # linear model
            VE_decays$est[VE_decays$model == "linear"] = 
              times_lin[3:6, -1, drop = FALSE] %*% coef(fit_lin)[-1]
            VE_decays$var[VE_decays$model == "linear"] = 
              apply(times_lin[3:6, -1, drop = FALSE], 1, 
                    function(x) t(x) %*% vcov(fit_lin)[-1, -1, drop = F] %*% x)
            
            # quadratic model
            VE_decays$est[VE_decays$model == "quadratic"] = 
              times_quad[3:6, -1, drop = FALSE] %*% coef(fit_quad)[-1]
            VE_decays$var[VE_decays$model == "quadratic"] = 
              apply(times_quad[3:6, -1, drop = FALSE], 1, 
                    function(x) t(x) %*% vcov(fit_quad)[-1, -1] %*% x)
            
            # pspline model
            VE_decays$est[VE_decays$model == "pspline"] = 
              apply(times_ps[3:6, -1, drop = F], 1, 
                    function(x) (x - centvec[-1]) %*% coef(fit_ps)[-1])
            
            VE_decays$var[VE_decays$model == "pspline"] = 
              apply(times_ps[3:6, -1, drop = F], 1, 
                    function(x) t(x - centvec[-1]) %*% vcov(fit_ps)[-1,-1] %*% (x - centvec[-1]))
            
            # return VE_ests
            return(list(par_ests = par_ests, 
                        VE_ests = VE_ests, 
                        VE_decays = VE_decays,
                        n_cross = n_cross,
                        t_cross = t_cross,
                        event_counts = event_counts))
          }

stopImplicitCluster()

saveRDS(trials, file = paste0("pure_sim_", simnum,".Rds"))

closeCluster(cl)
mpi.quit()
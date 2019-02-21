library(trophToxNTO)
## --------------------------- TIME ------------------------
times <- seq(from=0, to = 50*365, length = 1e4) ##  30*5 = 5 months
## --------------------------- EVENEMENTS -------------------
root.event=condition(Vole.Level=500,Broma.Level=20*50)[[1]]
event.func=condition(Vole.Level=500,Broma.Level=20*50)[[2]]
## ------------------------- RUN
params = pars(n=1)
init.p <- init(params$n)
# # Set fix parameters
# range_param = rbind(unlist(fix_pVar_pars(params,-0.1)), unlist(fix_pVar_pars(params,0.1)))
# ls_params = as.list(as.data.frame(range_param))
# ls_params_expand <- expand.grid(ls_params[2:ncol(range_param)])

# # Random parameters
ls_params_expand <- lapply(1:1e3, function(i) rand_pVar_pars(params))

ls_randSens <- lapply(1:length(ls_params_expand), function(i){
  tryCatch({
    model_output = ode(y=init.p,
                       times,
                       ode.broma.tritrophic,
                       parms = ls_params_expand[[i]],
                       events = list(func = event.func, root = TRUE),
                       rootfun = root.event)
    df_temp <- tab_sensitivity(model_output, ls_params_expand[[i]])
    return(df_temp)},
    error=function(e) NULL)
})
# mat_randSens <- do.call("rbind", ls_randSens)
df_randSens <- dplyr::bind_rows(ls_randSens)
# save(df_randSens, file = paste0("data/df_randSens_", Sys.time(), "_", round(runif(1,1e10,1e11)), ".rda"))
save(df_randSens, file = paste0("df_randSens_", Sys.time(), "_", round(runif(1,1e10,1e11)), ".rda"))

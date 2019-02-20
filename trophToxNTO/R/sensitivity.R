#' @title Table of sensitivity for set, min, max values
#'
#' @descritpiont Take the output of the ode integration with deSolve and return a data.frame
#' with a set of measure on
#' @param model_output the output of the ode solver
#' @param params Parameter of the system
#' @param key specific name of the parameters to consider
#'
#' @return A list of data.frame
#'
#' @importFrom zoo rollmean
#'
#' @export
#'
#'
tab_sensitivity <- function(model_output, params, tab_ref){

  tab_set <- table_df(model_output)
  df_params <- as.data.frame(params)
  df_params$Mean_V_set <- mean(tab_set$D.V)

  df_params$Mean_M_set <- mean(tab_set$D.M)

  df_params$Mean_F_set <- mean(tab_set$D.F)

  df_params$V_inf50_set <- nrow(tab_set[tab_set$D.V <50,]) / nrow(tab_set)

  df_params$nrb_AR_set <- nbr_spike(tab_set$C.soil)

  tab_mortality_set <- Nbr.mortality(tab_set,params,group.num=1)

  df_params$M_50AR_set <- sum(tab_mortality_set$Percent.M.mort.B > 0.5) / nrow(tab_mortality_set)

  df_params$M_50natural_set <- sum(tab_mortality_set$Percent.M.mort.Nm > 0.5) / nrow(tab_mortality_set)

  return(df_params)
}



tab_SSE <- function(model_output, params, tab_ref){
  tab_set <- table_df(model_output)
  mat_diff <- tab_ref[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")] - tab_set[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")]
  params$SSE <- sum(mat_diff^2)/(dim(mat_diff)[1] * dim(mat_diff)[2])
  return(params)
}

#' @title Table of sensitivity for set, min, max values
#'
#' @descritpiont Take the output of the ode integration with deSolve and return a data.frame
#' with a set of measure on
#' @param model_output the output of the ode solver
#' @param params Parameter of the system
#' @param key specific name of the parameters to consider
#'
#' @return A list of data.frame
#'
#' @importFrom zoo rollmean
#'
#' @export
#'
#'
tab_sensitivity2 <- function(model_output, params, key = "r.V"){

  df_out = data.frame(parameter = key)

  df_out$param_min = min(params[[key]])
  df_out$param_max = min(params[[key]])

  tab <- table_df(model_output)

  tab_set <- tab[tab$group == 1,]
  tab_min <- tab[tab$group == 2,]
  tab_max <- tab[tab$group == 3,]

  mat_min_diff <- tab_set[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")] - tab_min[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")]
  mat_max_diff <- tab_set[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")] - tab_max[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")]

  df_out$SSE_min <- sum(mat_min_diff^2)/(dim(mat_min_diff)[1] * dim(mat_min_diff)[2])
  df_out$SSE_max <- sum(mat_max_diff^2)/(dim(mat_max_diff)[1] * dim(mat_max_diff)[2])

  df_out$Mean_V_set <- mean(tab_set$D.V)
  df_out$Mean_V_min <- mean(tab_min$D.V)
  df_out$Mean_V_max <- mean(tab_max$D.V)

  df_out$Mean_M_set <- mean(tab_set$D.V)
  df_out$Mean_M_min <- mean(tab_min$D.V)
  df_out$Mean_M_max <- mean(tab_max$D.V)

  df_out$Mean_F_set <- mean(tab_set$D.V)
  df_out$Mean_F_min <- mean(tab_min$D.V)
  df_out$Mean_F_max <- mean(tab_max$D.V)

  df_out$V_inf50_set <- nrow(tab_set[tab_set$D.V <50,]) / nrow(tab_set)
  df_out$V_inf50_min <- nrow(tab_min[tab_min$D.V <50,]) / nrow(tab_min)
  df_out$V_inf50_max <- nrow(tab_max[tab_max$D.V <50,]) / nrow(tab_max)

  df_out$nrb_AR_set <- nbr_spike(tab_set$C.soil)
  df_out$nrb_AR_min <- nbr_spike(tab_min$C.soil)
  df_out$nrb_AR_max <- nbr_spike(tab_max$C.soil)

  tab_mortality_set = Nbr.mortality(tab,params,group.num=1)
  tab_mortality_min = Nbr.mortality(tab,params,group.num=2)
  tab_mortality_max = Nbr.mortality(tab,params,group.num=3)

  df_out$M_50AR_set <- sum(tab_mortality_set$Percent.M.mort.B > 0.5) / nrow(tab_mortality_set)
  df_out$M_50AR_min <- sum(tab_mortality_min$Percent.M.mort.B > 0.5) / nrow(tab_mortality_min)
  df_out$M_50AR_max <- sum(tab_mortality_max$Percent.M.mort.B > 0.5) / nrow(tab_mortality_max)

  df_out$M_50natural_set <- sum(tab_mortality_set$Percent.M.mort.Nm > 0.5) / nrow(tab_mortality_set)
  df_out$M_50natural_min <- sum(tab_mortality_min$Percent.M.mort.Nm > 0.5) / nrow(tab_mortality_min)
  df_out$M_50natural_max <- sum(tab_mortality_max$Percent.M.mort.Nm > 0.5) / nrow(tab_mortality_max)

  return(df_out)
}


#' @title random set of parameters
#'
#' @description random set of parameters
#'
#' @param params set of parameter as reference
#' @param p rate of variation min max in auniform distribution compared to reference
#'
#' @export
#'
rand_pVar_pars  <- function(params, p=0.1){
  # join parameters
  a.M = rand_pVar(params$a.M, p)
  a.Md = a.M
  a.F=rand_pVar(params$a.F, p)
  a.Fd=a.F
  return(list(
    n=1, # number of parallel simulation
    r.V = rand_pVar(params$r.V, p), # reproduction rate of voles - Nbr/days
    K.V = rand_pVar(params$K.V, p), # carrying capacity of voles - Nbr
    a.M = a.M,
    a.Md = a.M,
    a.F=a.F,
    a.Fd=a.F, # attack rate of Mustelid and Red fox - Nbr/day
    h.M=rand_pVar(params$h.M, p),
    h.F=rand_pVar(params$h.F, p), # handling time of Mustelid and Red fox - day
    m.M=rand_pVar(params$m.M, p),
    d=rand_pVar(params$d),
    #g.MF=g.MF, # consumption rate of Red fox over Mustelid
    a.FM=rand_pVar(params$a.FM, p), # attack rate of foxes on mustelids
    r.F=rand_pVar(params$r.F, p), # reproduction rate of Red fox - Nbr/days
    K.F=rand_pVar(params$K.F, p), # carrying capacity of Red fox - Nbr
    e.M=rand_pVar(params$e.M, p), # Energy conversion
    #
    LD.50.V=rand_pVar(params$LD.50.V, p),
    LD.50.M=rand_pVar(params$LD.50.M, p),
    LD.50.F=rand_pVar(params$LD.50.F, p),
    days.to.die.V=rand_pVar(params$days.to.die.V, p),
    days.to.die.M=rand_pVar(params$days.to.die.M, p),
    days.to.die.F=rand_pVar(params$days.to.die.F, p),
    H.V=rand_pVar(params$H.V, p),
    H.M=rand_pVar(params$H.M, p),
    H.F=rand_pVar(params$H.F, p),
    # Mean biomass of one individual populations - grams
    B.V = rand_pVar(params$B.V, p) , ## voles
    B.M = rand_pVar(params$B.M, p) , ## mustelid
    B.F = rand_pVar(params$B.F, p), ## fox
    max.intake.C=rand_pVar(params$max.intake.C, p), # daily intake: 6 ppm have been measured in vole after 4 à 6 jours : here 5 days
    D50.intake.C=rand_pVar(params$D50.intake.C, p), # at max.intake/2
    # Dynamic of contaminant VOIR Sage_2008 et Grolleau.
    eta.M = rand_pVar(params$eta.M, p),
    eta.F = rand_pVar(params$eta.F, p), ## Absorbtion rate
    k.out.V=rand_pVar(params$k.out.V, p),
    k.out.M=rand_pVar(params$k.out.M, p),
    k.out.F=rand_pVar(params$k.out.F, p), ## excretion rate
    k.0 = rand_pVar(params$k.0, p) # persitence broma in environment
  ))
}



#' @title fix set of parameters
#'
#' @description fix set of parameters
#'
#' @param params set of parameter as reference
#' @param p rate of variation compared to reference
#'
#' @export
#'
fix_pVar_pars  <- function(params, p){
  # join parameters
  a.M = fix_pVar(params$a.M, p)
  a.Md = a.M
  a.F=fix_pVar(params$a.F, p)
  a.Fd=a.F
  return(list(
    n=1, # number of parallel simulation
    r.V = fix_pVar(params$r.V, p), # reproduction rate of voles - Nbr/days
    K.V = fix_pVar(params$K.V, p), # carrying capacity of voles - Nbr
    a.M = a.M,
    a.Md = a.M,
    a.F=a.F,
    a.Fd=a.F, # attack rate of Mustelid and Red fox - Nbr/day
    h.M=fix_pVar(params$h.M, p),
    h.F=fix_pVar(params$h.F, p), # handling time of Mustelid and Red fox - day
    m.M=fix_pVar(params$m.M, p),
    d=fix_pVar(params$d, p),
    #g.MF=g.MF, # consumption rate of Red fox over Mustelid
    a.FM=fix_pVar(params$a.FM, p), # attack rate of foxes on mustelids
    r.F=fix_pVar(params$r.F, p), # reproduction rate of Red fox - Nbr/days
    K.F=fix_pVar(params$K.F, p), # carrying capacity of Red fox - Nbr
    e.M=fix_pVar(params$e.M, p), # Energy conversion
    #
    LD.50.V=fix_pVar(params$LD.50.V, p),
    LD.50.M=fix_pVar(params$LD.50.M, p),
    LD.50.F=fix_pVar(params$LD.50.F, p),
    days.to.die.V=fix_pVar(params$days.to.die.V, p),
    days.to.die.M=fix_pVar(params$days.to.die.M, p),
    days.to.die.F=fix_pVar(params$days.to.die.F, p),
    H.V=fix_pVar(params$H.V, p),
    H.M=fix_pVar(params$H.M, p),
    H.F=fix_pVar(params$H.F, p),
    # Mean biomass of one individual populations - grams
    B.V = fix_pVar(params$B.V, p) , ## voles
    B.M = fix_pVar(params$B.M, p) , ## mustelid
    B.F = fix_pVar(params$B.F, p), ## fox
    max.intake.C=fix_pVar(params$max.intake.C, p), # daily intake: 6 ppm have been measured in vole after 4 à 6 jours : here 5 days
    D50.intake.C=fix_pVar(params$D50.intake.C, p), # at max.intake/2
    # Dynamic of contaminant VOIR Sage_2008 et Grolleau.
    eta.M = fix_pVar(params$eta.M, p),
    eta.F = fix_pVar(params$eta.F, p), ## Absorbtion rate
    k.out.V=fix_pVar(params$k.out.V, p),
    k.out.M=fix_pVar(params$k.out.M, p),
    k.out.F=fix_pVar(params$k.out.F, p), ## excretion rate
    k.0 = fix_pVar(params$k.0, p) # persitence broma in environment
  ))
}




#' @title sequence of set, min, max parameter
#'
#' @description sequence of set, min, max parameter
#'
#' @export
#'
varPercent_seq <- function(x, p, n){
  if(p<0 || p>1) stop()
  return( c(x*(1-p), x, x * (1+p)) )
}

AUC <- function(x, y){
  sum(diff(x)*rollmean(y,2))
}

### INTERNAL
#' @export
nbr_spike <- function(y, k=10){
  y_smooth <- zoo::rollmean(y, k=k, na.pad = TRUE, align = "right") # 50*365 / 1e4 * k = nbr days
  # plot(x,y_smooth, type = "l")
  diff_y <- diff(y_smooth)
  diff_y_up1 <- c(diff_y[2:length(diff_y)], NA)
  bool_posTOneg <- (diff_y > 0 & diff_y_up1 <= 0)
  # points(x,c(bool_posTOneg*400, NA), col  = "red")
  return(sum(bool_posTOneg[!is.na(bool_posTOneg)]))
}

# Internal
rand_pVar <- function(x,p=0.1){
  # res <- runif(1, x*(1-p), x*(1+p))
  x_min = x*(1-p) ; x_max = x*(1+p)
  res <- rbeta(1,0.2,0.2)*(x_max - x_min)+x_min
  return(res)
}
fix_pVar <- function(x,p=-0.1){
  return(x*(1+p))
}

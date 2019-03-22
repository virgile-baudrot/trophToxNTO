#' @title Table of sensitivity for set, min, max values
#'
#' @description Take the output of the ode integration with deSolve and return a data.frame
#' with a set of measure on
#' @param model_output the output of the ode solver
#' @param parameters Parameter of the system
#' @param key specific name of the parameters to consider
#'
#' @return A list of data.frame
#'
#' @importFrom zoo rollmean
#'
#' @export
#'
#'
tab_sensitivity <- function(model_output, parameters){

  tab_set <- table_df(model_output)
  df_parameters <- as.data.frame(parameters)
  df_parameters$Mean_V <- mean(tab_set$D.V, na.rm = T)

  df_parameters$Mean_M <- mean(tab_set$D.M, na.rm = T)

  df_parameters$Mean_F <- mean(tab_set$D.F, na.rm = T)

  df_parameters$V_inf50 <- nrow(tab_set[tab_set$D.V <50,]) / nrow(tab_set)

  df_parameters$nrb_AR <- nbr_spike(tab_set$C.soil)

  tab_mortality <- Nbr.mortality(tab_set,parameters,group.num=1)

  df_parameters$M_50AR <- sum(tab_mortality$Percent.M.mort.B > 0.5) / nrow(tab_mortality)

  df_parameters$M_50natural <- sum(tab_mortality$Percent.M.mort.Nm > 0.5) / nrow(tab_mortality)

  return(df_parameters)
}



tab_SSE <- function(model_output, parameters, tab_ref){
  tab_set <- table_df(model_output)
  mat_diff <- tab_ref[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")] - tab_set[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")]
  parameters$SSE <- sum(mat_diff^2)/(dim(mat_diff)[1] * dim(mat_diff)[2])
  return(parameters)
}

#' @title Table of sensitivity for set, min, max values
#'
#' @description Take the output of the ode integration with deSolve and return a data.frame
#' with a set of measure on
#' @param model_output the output of the ode solver
#' @param parameters Parameter of the system
#' @param key specific name of the parameters to consider
#'
#' @return A list of data.frame
#'
#' @importFrom zoo rollmean
#'
#' @export
#'
#'
tab_sensitivity2 <- function(model_output, parameters, key = "r.V"){

  df_out = data.frame(parameter = rep(key, 3),
                      range = c("set", "min", "max"))

  df_out$param = c(NA, min(parameters[[key]]), max(parameters[[key]]))

  tab <- table_df(model_output)

  tab_set <- tab[tab$group == 1,]
  tab_min <- tab[tab$group == 2,]
  tab_max <- tab[tab$group == 3,]

  mat_min_diff <- tab_set[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")] - tab_min[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")]
  mat_max_diff <- tab_set[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")] - tab_max[, c("C.F", "C.M", "C.V", "D.F", "D.M", "D.V", "D.Vd")]

  df_out$SSE <- c(NA,
                  sum(mat_min_diff^2)/(dim(mat_min_diff)[1] * dim(mat_min_diff)[2]),
                  sum(mat_max_diff^2)/(dim(mat_max_diff)[1] * dim(mat_max_diff)[2]))

  df_out$Mean_V <- c(mean(tab_set$D.V),
                     mean(tab_min$D.V),
                     mean(tab_max$D.V))

  df_out$Mean_M <- c(mean(tab_set$D.M),
                     mean(tab_min$D.M),
                     mean(tab_max$D.M))

  df_out$Mean_F <- c(mean(tab_set$D.F),
                     mean(tab_min$D.F),
                     mean(tab_max$D.F))

  df_out$V_inf50 <- c( nrow(tab_set[tab_set$D.V <50,]) / nrow(tab_set),
                       nrow(tab_min[tab_min$D.V <50,]) / nrow(tab_min),
                       nrow(tab_max[tab_max$D.V <50,]) / nrow(tab_max))

  df_out$nrb_AR <- c(nbr_spike(tab_set$C.soil),
                     nbr_spike(tab_min$C.soil),
                     nbr_spike(tab_max$C.soil))

  tab_mortality_set = Nbr.mortality(tab,parameters,group.num=1)
  tab_mortality_min = Nbr.mortality(tab,parameters,group.num=2)
  tab_mortality_max = Nbr.mortality(tab,parameters,group.num=3)

  df_out$M_50AR <- c(sum(tab_mortality_set$Percent.M.mort.B > 0.5) / nrow(tab_mortality_set),
                     sum(tab_mortality_min$Percent.M.mort.B > 0.5) / nrow(tab_mortality_min),
                     sum(tab_mortality_max$Percent.M.mort.B > 0.5) / nrow(tab_mortality_max) )

  df_out$M_50natural <- c(sum(tab_mortality_set$Percent.M.mort.Nm > 0.5) / nrow(tab_mortality_set),
                          sum(tab_mortality_min$Percent.M.mort.Nm > 0.5) / nrow(tab_mortality_min),
                          sum(tab_mortality_max$Percent.M.mort.Nm > 0.5) / nrow(tab_mortality_max) )

  return(df_out)
}


#' @title random set of parameters
#'
#' @description random set of parameters
#'
#' @param parameters set of parameter as reference
#' @param p rate of variation min max in auniform distribution compared to reference
#'
#' @export
#'
rand_pVar_pars  <- function(parameters, p=0.1){

  if(length(p) == 1) p = rep(p, 32)
  # join parameters
  a.M = rand_pVar(parameters$a.M, p[3])
  a.Md = a.M
  a.F=rand_pVar(parameters$a.F, p[4])
  a.Fd=a.F
  return(list(
    n=1, # number of parallel simulation
    r.V = rand_pVar(parameters$r.V, p[1]), # reproduction rate of voles - Nbr/days
    K.V = rand_pVar(parameters$K.V, p[2]), # carrying capacity of voles - Nbr
    a.M = a.M,
    a.Md = a.M,
    a.F=a.F,
    a.Fd=a.F, # attack rate of Mustelid and Red fox - Nbr/day
    h.M=rand_pVar(parameters$h.M, p[5]),
    h.F=rand_pVar(parameters$h.F, p[6]), # handling time of Mustelid and Red fox - day
    m.M=rand_pVar(parameters$m.M, p[7]),
    d=rand_pVar(parameters$d, p[8]),
    #g.MF=g.MF, # consumption rate of Red fox over Mustelid
    a.FM=rand_pVar(parameters$a.FM, p[9]), # attack rate of foxes on mustelids
    r.F=rand_pVar(parameters$r.F, p[10]), # reproduction rate of Red fox - Nbr/days
    K.F=rand_pVar(parameters$K.F, p[11]), # carrying capacity of Red fox - Nbr
    e.M=rand_pVar(parameters$e.M, p[12]), # Energy conversion
    #
    LD.50.V=rand_pVar(parameters$LD.50.V, p[13]),
    LD.50.M=rand_pVar(parameters$LD.50.M, p[14]),
    LD.50.F=rand_pVar(parameters$LD.50.F, p[15]),
    days.to.die.V=rand_pVar(parameters$days.to.die.V, p[16]),
    days.to.die.M=rand_pVar(parameters$days.to.die.M, p[17]),
    days.to.die.F=rand_pVar(parameters$days.to.die.F, p[18]),
    H.V=rand_pVar(parameters$H.V, p[19]),
    H.M=rand_pVar(parameters$H.M, p[20]),
    H.F=rand_pVar(parameters$H.F, p[21]),
    # Mean biomass of one individual populations - grams
    B.V = rand_pVar(parameters$B.V, p[22]) , ## voles
    B.M = rand_pVar(parameters$B.M, p[23]) , ## mustelid
    B.F = rand_pVar(parameters$B.F, p[24]), ## fox
    max.intake.C=rand_pVar(parameters$max.intake.C, p[25]), # daily intake: 6 ppm have been measured in vole after 4 à 6 jours : here 5 days
    D50.intake.C=rand_pVar(parameters$D50.intake.C, p[26]), # at max.intake/2
    # Dynamic of contaminant VOIR Sage_2008 et Grolleau.
    eta.M = rand_pVar(parameters$eta.M, p[27]),
    eta.F = rand_pVar(parameters$eta.F, p[28]), ## Absorbtion rate
    k.out.V=rand_pVar(parameters$k.out.V, p[29]),
    k.out.M=rand_pVar(parameters$k.out.M, p[30]),
    k.out.F=rand_pVar(parameters$k.out.F, p[31]), ## excretion rate
    k.0 = rand_pVar(parameters$k.0, p[32]) # persitence broma in environment
  ))
}



#' @title fix set of parameters
#'
#' @description fix set of parameters
#'
#' @param parameters set of parameter as reference
#' @param p rate of variation compared to reference
#'
#' @export
#'
fix_pVar_pars  <- function(parameters, p){
  # join parameters
  a.M = fix_pVar(parameters$a.M, p)
  a.Md = a.M
  a.F=fix_pVar(parameters$a.F, p)
  a.Fd=a.F
  return(list(
    n=1, # number of parallel simulation
    r.V = fix_pVar(parameters$r.V, p), # reproduction rate of voles - Nbr/days
    K.V = fix_pVar(parameters$K.V, p), # carrying capacity of voles - Nbr
    a.M = a.M,
    a.Md = a.M,
    a.F=a.F,
    a.Fd=a.F, # attack rate of Mustelid and Red fox - Nbr/day
    h.M=fix_pVar(parameters$h.M, p),
    h.F=fix_pVar(parameters$h.F, p), # handling time of Mustelid and Red fox - day
    m.M=fix_pVar(parameters$m.M, p),
    d=fix_pVar(parameters$d, p),
    #g.MF=g.MF, # consumption rate of Red fox over Mustelid
    a.FM=fix_pVar(parameters$a.FM, p), # attack rate of foxes on mustelids
    r.F=fix_pVar(parameters$r.F, p), # reproduction rate of Red fox - Nbr/days
    K.F=fix_pVar(parameters$K.F, p), # carrying capacity of Red fox - Nbr
    e.M=fix_pVar(parameters$e.M, p), # Energy conversion
    #
    LD.50.V=fix_pVar(parameters$LD.50.V, p),
    LD.50.M=fix_pVar(parameters$LD.50.M, p),
    LD.50.F=fix_pVar(parameters$LD.50.F, p),
    days.to.die.V=fix_pVar(parameters$days.to.die.V, p),
    days.to.die.M=fix_pVar(parameters$days.to.die.M, p),
    days.to.die.F=fix_pVar(parameters$days.to.die.F, p),
    H.V=fix_pVar(parameters$H.V, p),
    H.M=fix_pVar(parameters$H.M, p),
    H.F=fix_pVar(parameters$H.F, p),
    # Mean biomass of one individual populations - grams
    B.V = fix_pVar(parameters$B.V, p) , ## voles
    B.M = fix_pVar(parameters$B.M, p) , ## mustelid
    B.F = fix_pVar(parameters$B.F, p), ## fox
    max.intake.C=fix_pVar(parameters$max.intake.C, p), # daily intake: 6 ppm have been measured in vole after 4 à 6 jours : here 5 days
    D50.intake.C=fix_pVar(parameters$D50.intake.C, p), # at max.intake/2
    # Dynamic of contaminant VOIR Sage_2008 et Grolleau.
    eta.M = fix_pVar(parameters$eta.M, p),
    eta.F = fix_pVar(parameters$eta.F, p), ## Absorbtion rate
    k.out.V=fix_pVar(parameters$k.out.V, p),
    k.out.M=fix_pVar(parameters$k.out.M, p),
    k.out.F=fix_pVar(parameters$k.out.F, p), ## excretion rate
    k.0 = fix_pVar(parameters$k.0, p) # persitence broma in environment
  ))
}




#' @title sequence of set, min, max parameter
#'
#' Sequence of set, min, max parameter
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
  res <- rbeta(1,0.01,0.01)*(x_max - x_min)+x_min
  return(res)
}
fix_pVar <- function(x,p=-0.1){
  return(x*(1+p))
}

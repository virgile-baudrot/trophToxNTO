#' @title General Solver for Ordinary Differential Equations
#'
#' @description function from package \link{deSolve}.
#' Solves a system of ordinary differential equations;
#'  a wrapper around the implemented ODE solvers
#'
#' @name ode
#' @rdname ode
#' @keywords internal
#' @export
#' @importFrom deSolve ode
NULL

#' @title Arrange multiple grobs on a page
#'
#' @description function from package \link{gridExtra}.
#' Set up a gtable layout to place multiple grobs on a page.
#'
#' @name grid.arrange
#' @rdname grid.arrange
#' @keywords internal
#' @export
#' @importFrom gridExtra grid.arrange
NULL

#' Dose-response curves for mortality
#'
#' Take a vecor containing : (i) the concentration in the body, (ii) the Lethal Dose of e species,
#'  (iii) the hill coeficient of the curve and (iv) the number of day to die
#'
#' @param C A numeric vector of the concentration in the individual
#' @param LD.50 A numeric value of the Lethal Dose 50 (50% of the individual die)
#' @param H.coef A numeric value coresponding to the Hill coefficient of the dose-response curve
#' @param days.death The number of days to die according to the LD50
#'
#'
#' @export
#'
#'
dose.resp.mort=function(C,LD.50,H.coef,days.death){
  return(
    1/days.death*(1-1/(1+(C/LD.50)^H.coef))
  )
}


#' Intake rate of contaminant for voles
#'
#' Take a vector containing
#'
#' @param C A numeric value/vector of contaminant concentration available
#' @param max.intake.C A numeric value of the maximal intake concentration
#' @param D50.intake A numeric value of the concentration in environment at 50% of max intake
#'
#' @return The Intake rate
#'
#' @export
#'
#'
intake.V=function(C, max.intake.C, D50.intake.C){
  return(max.intake.C*C/(D50.intake.C+C))
}

#' The system of differential equations to integrate with deSolve package
#'
#' Take Time, State, Pars as input, return derivative list
#'
#' @param Time
#' @param State
#' @param Pars
#'
#' @return A derivative list to integrate with deSolve package
#'
#' @export
#'
#'
ode.broma.tritrophic_2nd = function(Time, State, Pars) { # The model is a function
  with(as.list(c(State, Pars)), { # You need State and Pars in the next lines

    C.soil = State[1:n]
    C.V = State[(n+1):(2*n)]
    C.M = State[(2*n+1):(3*n)]
    C.F = State[(3*n+1):(4*n)]
    D.V = State[(4*n+1):(5*n)]
    D.Vd = State[(5*n+1):(6*n)]
    D.M = State[(6*n+1):(7*n)]
    D.F = State[(7*n+1):(8*n)]
    Cing.V = State[(8*n+1):(9*n)]
    Cing.M = State[(9*n+1):(10*n)]
    Cing.F = State[(10*n+1):(11*n)]

    # Contaminant concentration
    #C.soil=approx.broma.farmer.input(Time) # See afterward, for external events, we have to approximate for ode integration

    dC.soil = -k.0*C.soil - intake.V(C.soil,max.intake.C,D50.intake.C) * D.V * B.V  #because we work in kg of species, and parameterization is in g
    # dC.soil =  -k.0*C.soil

    dC.V = intake.V(C.soil,max.intake.C,D50.intake.C) - C.V*(k.out.V+r.V*(1-D.V/K.V))

    dC.M = eta.M * B.V/B.M * (a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*C.V -
      C.M*(k.out.M + e.M*B.V*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd)))

    dC.F = eta.F * (a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*
      ((a.F*D.V+a.Fd*D.Vd)*B.V/B.F * C.V+a.FM*D.M*B.M/B.F *C.M)/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M) - C.F*(k.out.F + r.F*(1-D.F/K.F))

    ### vole dynamics
    dV = r.V*D.V*(1-D.V/K.V) -
      dose.resp.mort(Cing.V, LD.50.V,H.V,days.to.die.V)*D.V - # broma mortality # 6 days to die
      a.M*D.V/(a.M*D.V+a.Md*D.Vd)*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*D.M - #mustelid functional response
      a.F*D.V/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F #fox functional response

    ### dead Vole dynamics:
    dVd = dose.resp.mort(Cing.V, LD.50.V,H.V,days.to.die.V)*D.V - d*D.Vd -# add vole dead by bromadiolone
      a.M*D.Vd/(a.M*D.V+a.Md*D.Vd)*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*D.M - #mustelid functional response
      a.F*D.Vd/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F

    ## mustelid dynamics
    dM = e.M*B.V*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*D.M- #numerical response
      m.M*D.M - # natural mortality
      dose.resp.mort(Cing.M,LD.50.M,H.M,days.to.die.M)*D.M - #broma poisoning
      a.FM*D.M/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F# death du to competition with foxes

    # fox dynamics
    dF = r.F*D.F*( 1 - D.F/K.F) - # logistic growth
      dose.resp.mort(Cing.F,LD.50.F,H.F,days.to.die.F)*D.F # broma poisoning

    # useless, but kee
    dCing.V =  max(Cing.M*( e.M*B.V*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd)) -  m.M -
      dose.resp.mort(Cing.M, LD.50.M, H.M, days.to.die.M) -
     a.FM*D.M / (a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/
      (1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F / D.M ), 0)
    # --- --- --- ---
    dCing.M =  B.V/B.M * (a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd)) * C.V -
      Cing.M*(1+e.M*B.V*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd)) -  m.M -
               dose.resp.mort(Cing.M, LD.50.M, H.M, days.to.die.M) -
               a.FM*D.M / (a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/
                 (1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F / D.M)

    dCing.F =  (a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*
       ((a.F*D.V+a.Fd*D.Vd)*B.V/B.F*C.V + a.FM*D.M*B.M/B.F*C.M)/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)  -
      Cing.F*(1+r.F*(1-D.F/K.F) - dose.resp.mort(Cing.F,LD.50.F,H.F,days.to.die.F))


    return(list(c(dC.soil,dC.V, dC.M,dC.F, dV, dVd, dM, dF,  dCing.V, dCing.M, dCing.F))) # return
  })
}

#' The system of differential equations to integrate with deSolve package
#'
#' Take Time, State, Pars as input, return derivative list
#'
#' @param Time
#' @param State
#' @param Pars
#'
#' @return A derivative list to integrate with deSolve package
#'
#' @export
#'
ode.broma.tritrophic = function(Time, State, Pars) { # The model is a function
  with(as.list(c(State, Pars)), { # You need State and Pars in the next lines

    C.soil = State[1:n]
    C.V = State[(n+1):(2*n)]
    C.M = State[(2*n+1):(3*n)]
    C.F = State[(3*n+1):(4*n)]
    D.V = State[(4*n+1):(5*n)]
    D.Vd = State[(5*n+1):(6*n)]
    D.M = State[(6*n+1):(7*n)]
    D.F = State[(7*n+1):(8*n)]
    Cing.V = State[(8*n+1):(9*n)]
    Cing.M = State[(9*n+1):(10*n)]
    Cing.F = State[(10*n+1):(11*n)]
    # Contaminant concentration
    #C.soil=approx.broma.farmer.input(Time) # See afterward, for external events, we have to approximate for ode integration

    dC.soil = -k.0*C.soil - intake.V(C.soil,max.intake.C,D50.intake.C) *D.V*B.V  #because we work in kg of species, and parameterization is in g #Sage et al. 2008 2007
    #dC.soil = 0

    dC.V = intake.V(C.soil,max.intake.C,D50.intake.C) - C.V*(k.out.V+r.V*(1-D.V/K.V))

    dC.M = eta.M * B.V/B.M * (a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*C.V -
      C.M*(k.out.M + e.M*B.V*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd)))

    dC.F = eta.F * (a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*
      ((a.F*D.V+a.Fd*D.Vd)*B.V/B.F * C.V+a.FM*D.M*B.M/B.F *C.M)/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M) - C.F*(k.out.F + r.F*(1-D.F/K.F))

    ### vole dynamics
    dV = r.V*D.V*(1-D.V/K.V) -
      dose.resp.mort(C.V,LD.50.V,H.V,days.to.die.V)*D.V - # broma mortality # 6 days to die
      a.M*D.V/(a.M*D.V+a.Md*D.Vd)*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*D.M - #mustelid functional response
      a.F*D.V/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F #fox functional response

    ### dead Vole dynamics:
    dVd = dose.resp.mort(C.V,LD.50.V,H.V,days.to.die.V)*D.V - d*D.Vd -# add vole dead by bromadiolone
      a.M*D.Vd/(a.M*D.V+a.Md*D.Vd)*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*D.M - #mustelid functional response
      a.F*D.Vd/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F

    ## mustelid dynamics
    dM = e.M*B.V*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*D.M- #numerical response
      m.M*D.M - # natural mortality
      dose.resp.mort(Cing.M,LD.50.M,H.M,days.to.die.M)*D.M - #broma poisoning
      a.FM*D.M/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F# death du to competition with foxes

    # fox dynamics
    dF = r.F*D.F*( 1 - D.F/K.F) - # logistic growth
      dose.resp.mort(Cing.F,LD.50.F,H.F,days.to.die.F)*D.F # broma poisoning

    dCing.V = intake.V(C.soil,max.intake.C,D50.intake.C) - C.V * (k.out.V+r.V*(1-D.V/K.V))
    dCing.M = eta.M * B.V/B.M * (a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*C.V -
      C.M*(k.out.M + e.M*B.V*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd)))
    dCing.F = eta.F * (a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2) *
      ((a.F*D.V+a.Fd*D.Vd)*B.V/B.F * C.V+a.FM*D.M*B.M/B.F *C.M)/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M) - C.F*(k.out.F + r.F*(1-D.F/K.F))

    return(list(c(dC.soil,dC.V, dC.M,dC.F, dV, dVd, dM, dF,dCing.V, dCing.M, dCing.F))) # return
  })
}

#' Set parameters
#'
#' Set parameters
#'
#'
#' @export
#'
#'
pars  <- function(
  n=2, # number of parallel simulations

  # maximal reproduction rate of voles - Nbr/days
  r.V =  log(2*600)/365, # population can grow from 1 to 600 over a year - see ???
  K.V = 600, # carrying capacity of voles - Nbr/hectars

  # attack rate of Mustelid and Red fox - Nbr/day
  a.M=1/30, a.Md=1/30, ## 1/2 selon Jones et al., 2011
  a.F=1/10, a.Fd=1/10,
  # handling time of Mustelid and Red fox - day
  h.M = 1/3.5, #1/nbr vole dayly = 50gr daily ~ 1 vole daily
  h.F = 1/6, # 300gr daily ~ 5-6 voles daily

  # mortality rate of mustelid
  m.M = 1/(0.8*365), # life expectancy of a mustelid (max is 6-7years, but 1.5-2 in average)
  d = 1/7, # 1/time of decomposition of vole in the field

  #g.MF=0.00005, # consumption/competition rate of Red fox over Mustelid
  a.FM = 1/10 * 80/300 , # attack rate of foxes on mustelids : a.F * BV/BM
  r.F = log(3)/365, # reproduction rate of Red fox - Nbr/days
  K.F = 0.03 , # carrying capacity of Red fox - Nbr/hectares 2-3 per km^2 = 0.02-0.03 /hectares : Ruette et al., 2003
  e.M = 0.03125, # 0.004 Energy conversion utilisé par Wang, Nagy, Gilg, Kuang, 2009
  # Mean biomass of one individual populations - grams
  B.V = 0.08 , ## voles
  B.M = 0.3 , ## mustelid : 200-600gr weasel ;
  B.F = 5.8, ## fox
  #
  max.intake.C = 6*1/5, # daily intake: 6 ppm have been measured in vole after 4 à 6 jours : here 5 days ; Sage et al. 2008
  D50.intake.C = 100, # at max.intake/2

  # LD.50.V=2*80*1e-3, LD.50.M=9.8*300*1e-3,, LD.50.F=0.15*5800*1e-3, # Thèse de Manon est de Michael Sage mg administré par Kg de masse corporelle
  # LD.50 - NOUS AVONS CONSTRUIT LE MODELE AVEC LES CONCENTRATIONS, pas avec les quantités totales !
  LD.50.V = 2, days.to.die.V = 6, # on est plutôt sur du 2-3 voir 6 ppm !!!
  LD.50.M = 2.4, days.to.die.M = 5,
  LD.50.F = 0.5, days.to.die.F = 5, # LD.50.F or 0.15 et 5

  H.V=4, H.M=4, H.F=2,
  # Dynamic of contaminant
  eta.M = 0.5,
  eta.F = 0.5, ## Absorbtion rate de 69% à 100% de 4 à 24h (voir thèse Manon Jacquot)
  k.out.V = 0.4, # voir Sage et al. 2008
  k.out.M = 0.6,
  k.out.F = 0.6, ## excretion rate
  k.0 = (0.106+0.057)/2 #Sage 2008 et 2007 : persistence of the molecule in the environment
  #C.soil= # concentration of contaminant in the soil
){
  return(list(
    n=n, # number of parallel simulation

    r.V=r.V, # reproduction rate of voles - Nbr/days
    K.V=K.V, # carrying capacity of voles - Nbr
    a.M=a.M, a.Md=a.Md, a.F=a.F, a.Fd=a.Fd, # attack rate of Mustelid and Red fox - Nbr/day
    h.M=h.M, h.F=h.F, # handling time of Mustelid and Red fox - day
    m.M=m.M,
    d=d,

    #g.MF=g.MF, # consumption rate of Red fox over Mustelid
    a.FM=a.FM, # attack rate of foxes on mustelids

    r.F=r.F, # reproduction rate of Red fox - Nbr/days
    K.F=K.F, # carrying capacity of Red fox - Nbr
    e.M=e.M, # Energy conversion
    #
    LD.50.V=LD.50.V, LD.50.M=LD.50.M, LD.50.F=LD.50.F,
    days.to.die.V=days.to.die.V,
    days.to.die.M=days.to.die.M,
    days.to.die.F=days.to.die.F,
    H.V=H.V,H.M=H.M,H.F=H.F,
    # Mean biomass of one individual populations - grams
    B.V = B.V , ## voles
    B.M = B.M , ## mustelid
    B.F = B.F, ## fox
    max.intake.C=max.intake.C, # daily intake: 6 ppm have been measured in vole after 4 à 6 jours : here 5 days
    D50.intake.C=D50.intake.C, # at max.intake/2
    # Dynamic of contaminant VOIR Sage_2008 et Grolleau.
    eta.M = eta.M, eta.F = eta.F, ## Absorbtion rate
    k.out.V=k.out.V, k.out.M=k.out.M, k.out.F=k.out.F, ## excretion rate
    k.0 = k.0 # persitence broma in environment
  ))
}

#' Parameters for the system of differential equations to integrate with deSolve package
#'
#' Take all parameters a vector and return a list
#'
#' @param n.simu number of simulation to take into account
#' @param C.soil Concentration broma available in the soil
#' @param C.V  Concentration broma in individual vole
#' @param C.M  Concentration broma in individual mustelids
#' @param C.F  Concentration broma in individual fox
#' @param D.V  Density of the vole population
#' @param D.Vd  Density of the dead vole sub-population
#' @param D.M  Density of the mustelid population
#' @param D.F  Density of the fox population
#'
#' @return A list of numeric containing the initial condition for the differential equations to integrate with deSolve package
#'
#' @export
#'
#'
#'
init  <- function(n.simu,C.soil=0, C.V=0, C.M=0, C.F=0,  D.V=100, D.Vd=0, D.M=0.05, D.F=0.025, Cing.V=0, Cing.M=0, Cing.F=0){
  n=n.simu
  if(length(C.soil)==1){C.soil_ = rep(C.soil,n)} else {C.soil_ = C.soil}
  if(length(C.V)==1){C.V_=rep(C.V,n)} else {C.V_ = C.V}
  if(length(C.M)==1){C.M_=rep(C.M,n)} else {C.M_ = C.M}
  if(length(C.F)==1){C.F_=rep(C.F,n)} else {C.F_ = C.F}
  if(length(D.V)==1){D.V_ =rep(D.V,n)} else {D.V_ = D.V}
  if(length(D.Vd)==1){D.Vd_=rep(D.Vd,n)} else {D.Vd_ = D.Vd}
  if(length(D.M)==1){D.M_=rep(D.M,n)} else {D.M_ = D.M}
  if(length(D.F)==1){D.F_ = rep(D.F,n) } else {D.F_ = D.F}
  if(length(Cing.V)==1){Cing.V_ = rep(Cing.V,n)} else {Cing.V_ = Cing.V}
  if(length(Cing.M)==1){Cing.M_ = rep(Cing.M,n)} else {Cing.M_ = Cing.M}
  if(length(Cing.F)==1){Cing.F_ = rep(Cing.F,n) } else {Cing.F_ = Cing.F}
  return(c(C.soil_ = C.soil_ , C.V_=C.V_, C.M_=C.M_, C.F_=C.F_,
           D.V_ =D.V_, D.Vd_=D.Vd_,  D.M_=D.M_, D.F_ = D.F_ , Cing.V_ = Cing.V_, Cing.M_ = Cing.M_, Cing.F_ =Cing.F_))
}




#' Condition (density of voles) for event to happen (event root in deSolve parckage) And the events (input of broma in the field)
#'
#' Take a numeric and return the trigger condition of the event
#' @param Vole.Level A numeric value defining the level of vole for treatment
#' @param Broma.Level A numeric value defining the quantity of broma at treatment
#'
#' @return A list of two elements : the trigger moment of the event, and the event
#'
#' @export
#'
#'
condition = function(Vole.Level,Broma.Level){
  root.event = function(Time, State, Pars){
    with(as.list(c(State, Pars)), { D.V = State[(4*n+1):(5*n)]
    round((D.V - Vole.Level),5) # V-50 = 0 that is when V=50
    })
  }

  event.func = function(Time, State, Pars){
    with(as.list(c(State, Pars)), {
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

      C.soil=(round((D.V - Vole.Level),5)==0)*Broma.Level #rep(Broma.Level,n) # ppm * kg

      C.V=C.V ;  C.M=C.M ; C.F=C.F ;  D.V = D.V ;  D.Vd= D.Vd ;  D.M= D.M ;  D.F = D.F ; Cing.V = Cing.V ; Cing.M = Cing.M ; Cing.F = Cing.F
      return(c(C.soil,C.V,C.M,C.F,D.V,D.Vd,D.M,D.F,Cing.V, Cing.M,Cing.F))
    })
  }

  return(list(root.event,event.func))
}

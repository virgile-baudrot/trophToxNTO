## ---- eval=FALSE---------------------------------------------------------
#  install.packages('devtools') # if not already installed
#  library('devtools') # load the library
#  github_install("vbaudrot/trophToxNTO/trophToxNTO") # install from 'master'

## ----loadLibrary---------------------------------------------------------
library(trophToxNTO)

## ----warning=FALSE-------------------------------------------------------
vole.level.test = 250
broma.level.test = 20*50

# Run condition function for scenario
scenario = condition(Vole.Level=vole.level.test,Broma.Level=broma.level.test)

root.event=scenario[[1]]
event.func=scenario[[2]]

## ----warning=FALSE-------------------------------------------------------
# Load package deSolve
# library(deSolve)

# number of parallel simulation
params=pars(n=1)

# times vector of integration
times <- seq(from=0, to=15*365, length.out = 1e3) #  30*5 = 5 months 

# Run the model

model.output    <- ode(times = times,
                       y = init(params$n),
                       func = ode.broma.tritrophic,
                       parms  = params,
                       events = list(func = event.func, root = TRUE),
                       rootfun = root.event)


## ----warning=FALSE,message=FALSE-----------------------------------------
graphique(params$n,times,model.output) 

## ----warning=FALSE-------------------------------------------------------
# number of parallel simulation
params=pars(n=2)

# Run the model
model.output    <- ode(times = times,
                       y = init(params$n),
                       func = ode.broma.tritrophic,
                       parms = pars(
                         n=params$n, # Do not foret to specify this !!!
                         K.V=seq(200,800,length=params$n)
                       ),
                       events = list(func = event.func, root = TRUE),
                       rootfun = root.event)

# Graphic output
graphique(params$n,times,model.output) 

## ----warning=FALSE-------------------------------------------------------
# number of parallel simulation
params=pars(n=2)

# Run the model
model.output    <- ode(times = times,
                       y = init(params$n,
                                D.F = c(0,0.2)
                                ),
                       func = ode.broma.tritrophic,
                       parms = params,
                       events = list(func = event.func, root = TRUE),
                       rootfun = root.event)

# Graphic output
graphique(params$n,times,model.output) 


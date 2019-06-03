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

## ----warning=FALSE,message=FALSE-----------------------------------------
graphique(params$n,times,model.output) 


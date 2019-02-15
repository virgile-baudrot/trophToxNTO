#' Check if R package is installed then load library
#'
#' Take a package x check if installed, install it if necessary and load it: Check if R package is installed then load library
#' @param x A character corresponding to the package to load
#'
#'
#' @export
#'
#'

pkgLoad <- function(x){
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
  #now load library
  library(x,character.only=TRUE)
}

#' Table to use in graphics
#'
#' @param model.output The output of the ode integration with deSolve
#'
#' @import dplyr
#'
#' @return A data.frame
#'
#' @export
#'
#'
table_df = function(model.output){
  pkgLoad("dplyr")
  pkgLoad("tidyr")
  pkgLoad("stringr")

  df=model.output%>%
    as.data.frame()%>%
    gather(compartiment,value,-time) %>%
    separate(compartiment,c("compartiment","group"),sep="_") %>%
    mutate(group= ifelse(group=="","1",group))%>%
    spread(compartiment,value)
}

#' Table to use in graphics
#'
#' @param model.output The output of the ode integration with deSolve
#'
#' @return A data.frame
#'
#' @export
#'
#'
tableau = function(model.output){
  pkgLoad("dplyr")
  pkgLoad("tidyr")
  pkgLoad("stringr")

  df.1=table_df(model.output)

  df.g.2 = df.1 %>%
    gather(compartiment,value,-c(time,group))%>%
    # Add columns and information
    #+ mutate(category = ifelse(grepl('C', compartiment), 'concentration', compartiment))
    mutate(category = ifelse(grepl('C', compartiment), 'concentration', 'density'))%>% #-
    mutate(sub.category =
             ifelse(grepl('C', compartiment),
                    # concentration
                    ifelse(grepl('C.soil',compartiment) , 'C.soil','conc.species'),
                    # density
                    ifelse(grepl('^[^C]*[MF][^C]*$',compartiment), 'predator', 'vole')
             )
    )

  return(df.g.2)

}

#' Graphic for results of the system
#'
#' Take the output of the ode integration with deSolve and return graphics
#' @param n A numeric value of the number of simulation
#' @param times A numeric vector of the period of integration of the differential equation
#' @param model.output The output of the ode integration with deSolve
#'
#' @return A graphic
#'
#' @export
#'
#'
graphique = function(n,times,model.output){
  # ------------------- Package required --------------------------------
  pkgLoad("ggplot2")
  pkgLoad("dplyr")
  pkgLoad("tidyr")
  pkgLoad("stringr")

  # ----------- Manipulation des tableaux ----------------
  df.g.2 = tableau(model.output)

  n.treat= attributes(model.output)$nroot
  #tot.treat= n.treat*max(df.1$C.soil)

  # ------------------------------- Colors

  # Color palette
  col.V =  '#09708F'; col.Vd = '#8FBECD' ; col.M =  '#5E8F09'; col.Fx = '#FFA200'
  col.C.soil =  '#C15757' ; col.C.V = '#487F8F' ; col.C.M = '#A9CF68' ; col.C.F =  '#FFD180'

  #cbPalette <- c(col.V, col.Vd, col.M, col.Fx, col.C.soil, col.C.V, col.C.M, col.C.F)
  cbPalette <- rev(c(col.V, col.Vd, col.M, col.Fx, col.C.soil, col.C.V, col.C.M, col.C.F))

  # -------------------------------- levels
  #levels(factor(df.g.2$sub.category)) # de haut en bas
  df.g.2$sub.category <- factor(df.g.2$sub.category, levels = c('C.soil', 'conc.species', 'predator', 'vole'))
  #levels(factor(df.g.2$compartiment))
  #df.g.2$compartiment <- factor(df.g.2$compartiment , levels= c('V','Vd','M','Fx','C.soil', 'C.V', 'C.M', 'C.F' ))
  df.g.2$compartiment <- factor(df.g.2$compartiment , levels= rev(c('D.V','D.Vd','D.M','D.F','C.soil', 'C.V', 'C.M', 'C.F' )))

  plot.general=df.g.2 %>%
    ggplot()+theme_bw()+
    # ggtitle(paste("Number treatments", n.treat
    #               #, "\n", "Amount treatments", tot.treat, "mg/ha"
    #               ))+
    theme(legend.position="top")+
    scale_colour_manual(values=cbPalette)+
    scale_x_continuous(
      breaks = seq(0,max(times),length=11),
      labels=round(seq(0,max(times)/365,length=11),1),
      name="time in year")+
    scale_y_continuous(name="densities (ind/ha) and concentrations (ppm)")+
    #scale_y_sqrt(name="densisties and concentration")+
    #scale_y_continuous(limits = c(0, max(df$V)))+
    geom_line(aes(
      x=time,
      y=round(value,20),
      #group=interaction(group,compartiment),
      linetype=group,
      colour=compartiment),
      size=0.5,alpha=0.8)+
    #facet_grid(sub.category  ~ .,scales = "free")
    facet_grid(compartiment  ~ .,scales = "free")

  return(plot.general)
}

#' Graphic for results of the system
#'
#' Take the output of the ode integration with deSolve and return graphics
#' @param n A numeric value of the number of simulation
#' @param times A numeric vector of the period of integration of the differential equation
#' @param model.output The output of the ode integration with deSolve
#'
#' @return A graphic
#'
#' @export
#'


graphique_2 = function(n,times,model.output){
  # ------------------- Package required --------------------------------
  pkgLoad("ggplot2")
  pkgLoad("dplyr")
  pkgLoad("tidyr")
  pkgLoad("stringr")

  # ----------- Manipulation des tableaux ----------------
  df=as.data.frame(model.output) # data.frame of the ode

  df.1=df%>%
    gather(group,value,-time)%>%
    separate(group,c("species","group"),sep="_")%>%
    spread(species,value)%>%
    mutate(group.ren= ifelse(group=="1","No Fox", "With Fox") )

  n.treat= attributes(model.output)$nroot
  #tot.treat= n.treat*max(df.1$C.soil)

  df.g.2 = df.1 %>%
    gather(compartiment,value,-c(time,group,group.ren))%>%
    # Add columns and information
    #mutate(category = ifelse(grepl('C', compartiment), 'concentration', 'density'))
    mutate(category = ifelse(grepl('C', compartiment), 'concentration', compartiment))

  # Color palette
  col.V =  '#09708F'; col.Vd = '#8FBECD' ; col.M =  '#5E8F09'; col.Fx = '#FFA200'
  col.C.soil =  '#C15757' ; col.C.V = '#487F8F' ; col.C.M = '#A9CF68' ; col.C.F =  '#FFD180'


  #cbPalette <- c(col.V, col.Vd, col.M, col.Fx, col.C.soil, col.C.V, col.C.M, col.C.F)
  cbPalette <- rev(c(col.V, col.Vd, col.M, col.Fx, col.C.soil, col.C.V, col.C.M, col.C.F))

  # -------------------------------- levels
  #levels(factor(df.g.2$category)) # de haut en bas
  df.g.2$category <- factor(df.g.2$category, levels = c('concentration', 'D.F', 'D.M', 'D.Vd', 'D.V'))

  #levels(factor(df.g.2$sub.category)) # de haut en bas
  #df.g.2$sub.category <- factor(df.g.2$sub.category, levels = c('C.soil', 'conc.species', 'predator', 'vole'))
  #levels(factor(df.g.2$compartiment))
  df.g.2$compartiment <- factor(df.g.2$compartiment , levels= rev(c('D.V','D.Vd','D.M','D.F','C.soil', 'C.V', 'C.M', 'C.F' )))

  plot.general=df.g.2 %>%
    ggplot()+theme_bw()+
    ggtitle(paste("Number treatments", n.treat
                  #, "\n","Amount treatments", tot.treat, "mg/ha"
    ))+
    theme(legend.position="top")+
    scale_colour_manual(values=cbPalette)+
    scale_x_continuous(
      breaks = seq(0,max(times),length=11),
      labels=round(seq(0,max(times)/365,length=11),1),
      name="time in year")+
    scale_y_continuous(name="densities (ind/ha) and concentrations (ppm)")+
    geom_line(aes(
      x=time,
      y=round(value,20),
      #group=interaction(group,compartiment),
      linetype=group.ren,
      colour=compartiment),
      size=0.5,alpha=0.8)+
    facet_grid(category  ~ .,scales = "free")

  return(plot.general)
}



#' Graphic for results of the system
#'
#' Take the output of the ode integration with deSolve and return graphics
#' @param n A numeric value of the number of simulation
#' @param times A numeric vector of the period of integration of the differential equation
#' @param model.output The output of the ode integration with deSolve
#'
#' @return A graphic
#'
#' @export
#'


graphique_3 = function(n,times,model.output){
  # ------------------- Package required --------------------------------
  pkgLoad("ggplot2")
  pkgLoad("dplyr")
  pkgLoad("tidyr")
  pkgLoad("stringr")

  # ----------- Manipulation des tableaux ----------------
  df.1 <- table_df(model.output)

  n.treat= attributes(model.output)$nroot
  #tot.treat= n.treat*max(df.1$C.soil)

  df.g.2 = df.1 %>%
    gather(compartiment,value,-c(time,group))%>%
    # filter
    filter(compartiment %in% c("C.soil","D.V","D.M","D.F", "D.Vd"))

  # Color palette
  col.V =  '#09708F'; col.Vd =  '#000506'
  col.M =  '#5E8F09'; col.Fx = '#FFA200'
  col.C.soil =  '#C15757'

  cbPalette <- rev(c(col.C.soil, col.Vd, col.V, col.M, col.Fx))

  # -------------------------------- levels
  #levels(factor(df.g.2$compartiment)) # de haut en bas
  df.g.2$compartiment <- factor(df.g.2$compartiment, levels = c('D.F', 'D.M',  'D.V','D.Vd', 'C.soil'))

  plot.general=df.g.2 %>%
    ggplot()+theme_bw()+
    # ggtitle(paste("Number treatments", n.treat
    #               #, "\n","Amount treatments", tot.treat, "mg/ha"
    # ))+
    theme(legend.position="top")+
    scale_colour_manual(values=cbPalette)+
    scale_x_continuous(
      breaks = seq(0,max(times),length=11),
      labels=round(seq(0,max(times)/365,length=11),1),
      name="time in year")+
    scale_y_continuous(name="densities (ind/ha) and concentrations (ppm)")+
    geom_line(aes(
      x=time,
      y=round(value,20),
      #group=interaction(group,compartiment),
      linetype=group,
      colour=compartiment),
      size=0.5,alpha=0.8)+
    facet_grid(compartiment  ~ .,scales = "free")

  return(plot.general)
}

#' Repartition of the ingestion in mustelid and fox populations
#'
#' Take the output of the ode integration with deSolve and return a list of data.frame
#' @param table_df A data.frame of ode system
#' @param params Parameter of the system
#' @param group.num Group number to apply the computation
#'
#' @return A list of data.frame
#'
#' @export
#'
#'

Nbr.ingest=function(table_df,params,group.num=1){
  # -------------------- Munipulation Table ---------------
  df.mortality = dplyr::filter(table_df,group==group.num)

  with(c(params,df.mortality),{

    Nbr.V.ingest.M =  a.M * D.V/(a.M * D.V + a.Md * D.Vd) * (a.M * D.V + a.Md * D.Vd)/(1 + h.M * (a.M * D.V + a.Md * D.Vd)) * D.M

    Nbr.V.ingest.F = a.F * D.V/(a.F * D.V + a.Fd * D.Vd + a.FM * D.M) * (a.F * D.V + a.Fd * D.Vd + a.FM * D.M)^2/(1 + h.F * (a.F * D.V + a.Fd * D.Vd + a.FM * D.M)^2) * D.F

    Nbr.Vd.ingest.M = a.M * D.Vd/(a.M * D.V + a.Md * D.Vd) * (a.M * D.V + a.Md * D.Vd)/(1 + h.M * (a.M * D.V + a.Md * D.Vd)) * D.M

    Nbr.Vd.ingest.F = a.F * D.Vd/(a.F * D.V + a.Fd * D.Vd + a.FM * D.M) * (a.F * D.V + a.Fd * D.Vd + a.FM * D.M)^2/(1 + h.F * (a.F * D.V + a.Fd * D.Vd + a.FM * D.M)^2) * D.F

    Nbr.M.ingest.F = a.FM * D.M/(a.F * D.V + a.Fd * D.Vd + a.FM * D.M) * (a.F * D.V + a.Fd *  D.Vd + a.FM * D.M)^2/(1 + h.F * (a.F * D.V + a.Fd *D.Vd + a.FM * D.M)^2) * D.F

    sum.ingest.M = Nbr.V.ingest.M+Nbr.Vd.ingest.M
    sum.ingest.F=Nbr.V.ingest.F+Nbr.Vd.ingest.F+Nbr.M.ingest.F

    return(data.frame(time,
                      Nbr.V.ingest.M,
                      Nbr.Vd.ingest.M,
                      Nbr.V.ingest.F,
                      Nbr.Vd.ingest.F,
                      Nbr.M.ingest.F,
                      Percent.V.ingest.M  = Nbr.V.ingest.M  / sum.ingest.M,
                      Percent.Vd.ingest.M = Nbr.Vd.ingest.M / sum.ingest.M,
                      Percent.V.ingest.F  = Nbr.V.ingest.F  / sum.ingest.F,
                      Percent.Vd.ingest.F = Nbr.Vd.ingest.F / sum.ingest.F,
                      Percent.M.ingest.F  = Nbr.M.ingest.F  / sum.ingest.F
    ))
  })
}


#' Repartition of the mortality in vole and mustelid populations
#'
#' Take the output of the ode integration with deSolve and return a list of data.frame
#' @param table_df A data.frame of ode system
#' @param params Parameter of the system
#' @param group.num Group number to apply the computation
#'
#' @return A list of data.frame
#'
#' @export
#'
#'

Nbr.mortality=function(table_df,params,group.num=1){
  # -------------------- Munipulation Table ---------------
  df.mortality=dplyr::filter(table_df,group==group.num)

  with(c(params,df.mortality),{
    Nbr.V.mort.B =  dose.resp.mort(C.V,LD.50.V,H.V,days.to.die.V)*D.V
    Nbr.V.mort.MP = a.M*D.V/(a.M*D.V+a.Md*D.Vd)*(a.M*D.V+a.Md*D.Vd)/(1+h.M*(a.M*D.V+a.Md*D.Vd))*D.M
    Nbr.V.mort.FP =  a.F*D.V/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F
    Nbr.M.mort.B =  dose.resp.mort(C.M,LD.50.M,H.M,days.to.die.M)*D.M
    Nbr.M.mort.FP = a.FM*D.M/(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2/(1+h.F*(a.F*D.V+a.Fd*D.Vd+a.FM*D.M)^2)*D.F
    Nbr.M.mort.Nm = m.M*D.M

    sum.mort.V = Nbr.V.mort.B + Nbr.V.mort.MP +Nbr.V.mort.FP
    sum.mort.M =  Nbr.M.mort.B +  Nbr.M.mort.FP + Nbr.M.mort.Nm
    return(data.frame(time,
                      Nbr.V.mort.B,
                      Nbr.V.mort.FP,
                      Nbr.V.mort.MP,
                      Nbr.M.mort.B,
                      Nbr.M.mort.FP,
                      Nbr.M.mort.Nm,
                      Percent.V.mort.B = Nbr.V.mort.B/sum.mort.V,
                      Percent.V.mort.FP = Nbr.V.mort.FP/sum.mort.V,
                      Percent.V.mort.MP = Nbr.V.mort.MP/sum.mort.V,
                      Percent.M.mort.B = Nbr.M.mort.B/sum.mort.M,
                      Percent.M.mort.FP = Nbr.M.mort.FP/sum.mort.M,
                      Percent.M.mort.Nm = Nbr.M.mort.Nm/sum.mort.M
    ))
  })
}


#' Plot of the repartion of population ingested
#'
#' Take the output of the ode integration with deSolve and return a graphic
#' @param df.Nbr.ingest A data.frame of ode system
#'
#' @return A graphic
#'
#' @export
#'

plotIngest=function(df.Nbr.ingest){
  # -------------------- Manipulation Table ---------------
  df.use=df.Nbr.ingest %>%
    gather(group,value,-c(time))%>%
    mutate(category = ifelse(grepl('F', group), 'Fox', 'Mustelid'))%>%
    mutate(dimension = ifelse(grepl('Nbr', group), 'Number', 'Proportion'))

  df.use.M = df.use%>%
    dplyr::filter(category=="Mustelid")
  df.use.F = df.use%>%
    dplyr::filter(category=="Fox")

  #levels(factor( df.use.M$group))
  #levels(factor( df.use.F$group))
  # Color palette
  col.V =  '#09708F' ; col.Vd = '#8FBECD' ; col.M =  '#5E8F09'

  cbPalette.M <- c(col.Vd, col.V,  # Nbr
                   col.Vd, col.V ) # Percent

  cbPalette.F <- c(col.M, col.Vd, col.V, # Nbr
                   col.M, col.Vd, col.V) # Percent


  # -------------------- Plot ---------------
  plot.Mustelid=df.use.M%>%
    ggplot()+theme_bw()+
    ggtitle("Repartition ingestion for Mustelids")+
    scale_x_continuous(
      breaks = seq(0,max(times),length=11),
      labels=round(seq(0,max(times)/365,length=11),1),
      name="time in year")+
    theme(legend.position = "top")+
    scale_fill_manual(values=cbPalette.M)+
    geom_area(aes(x=time,y=value,fill=group))+
    facet_grid(dimension~.,scale="free")

  plot.Fox=df.use.F%>%
    ggplot()+theme_bw()+
    ggtitle("Repartition ingestion for Foxes")+
    scale_x_continuous(
      breaks = seq(0,max(times),length=11),
      labels=round(seq(0,max(times)/365,length=11),1),
      name="time in year")+
    scale_fill_manual(values=cbPalette.F)+
    theme(legend.position = "top")+
    geom_area(aes(x=time,y=value,fill=group))+
    facet_grid(dimension~.,scale="free")

  return(list(plot.Mustelid,plot.Fox))
}

#' Plot of the repartion of population mortality
#'
#' Take the output of the ode integration with deSolve and return a graphic
#' @param df.Nbr.ingest A data.frame of ode system
#'
#' @return A graphic
#'
#'
#' @export
#'

plotMortality=function(df.Nbr.mortality){
  # ------------------- Package required --------------------------------
  pkgLoad("dplyr")
  pkgLoad("tidyr")
  pkgLoad("ggplot2")

  # -------------------- Munipulation Table ---------------
  df.use=df.Nbr.mortality %>%
    gather(group,value,-c(time))%>%
    mutate(category = ifelse(grepl('V', group), 'Vole', 'Mustelid'))%>%
    mutate(dimension = ifelse(grepl('Nbr', group), 'Number', 'Proportion'))

  df.use.V = df.use%>%
    dplyr::filter(category=="Vole")
  df.use.M = df.use%>%
    dplyr::filter(category=="Mustelid")

  #levels(factor( df.use.V$group))
  #levels(factor( df.use.M$group))
  # Color palette
  col.B =  '#C15757' ; col.M =  '#5E8F09' ;  col.Fx = '#FFA200'

  cbPalette <- c(col.B, col.Fx, col.M, # Nbr.mort
                 col.B, col.Fx, col.M) # Percent


  # -------------------- Plot ---------------
  plot.Vole=df.use.V%>%
    ggplot()+theme_bw()+
    ggtitle("Repartition mortality in Vole population")+
    scale_x_continuous(
      breaks = seq(0,max(times),length=11),
      labels=round(seq(0,max(times)/365,length=11),1),
      name="time in year")+
    theme(legend.position = "top")+
    scale_fill_manual(values=cbPalette)+
    geom_area(aes(x=time,y=value,fill=group))+
    facet_grid(dimension~.,scale="free")

  plot.Mustelid=df.use.M%>%
    ggplot()+theme_bw()+
    ggtitle("Repartition mortality in Mustelid population")+
    scale_x_continuous(
      breaks = seq(0,max(times),length=11),
      labels=round(seq(0,max(times)/365,length=11),1),
      name="time in year")+
    scale_fill_manual(values=cbPalette)+
    theme(legend.position = "top")+
    geom_area(aes(x=time,y=value,fill=group))+
    facet_grid(dimension~.,scale="free")
  return(list(plot.Vole,plot.Mustelid))
}

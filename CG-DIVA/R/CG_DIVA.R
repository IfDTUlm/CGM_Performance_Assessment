# Continuous Glucose Deviation Interval and Variability Analysis (CG-DIVA)

# For documentation see
# https://github.com/IfDTUlm/CGM_Performance_Assessment

# Created by: Institut für Diabetes-Technology Forschungs- und Entwichlungsgesellschaft mbH an der Universität Ulm
# Contact: cgm_performance@idt-ulm.de

# This is a free software and comes with ABSOLUTELY NO WARRANTY

# Load packages 
library(data.table)       # setDT and mutate
library(pracma)           # linspace function
library(dplyr)            # %>% functionality
library(tidyr)            # nest functionality
library(purrr)            # map functionality
library(rsample)          # bootstrapping


# Plotting
plotting <- function(df,RES,version,ylims=c(-80,80),s_max=25,figsize=c(30,15))
{
  
  dat_int <- 0.9
  
  # Setup figure
  # windows(width=figsize[1]/2.5,height=figsize[2]/2.5)
  
  # Deviation Intervals
  # -------------------------------------------------------------------------
  
  # dev.new()
  # Setup axis
  par(mfrow = c(1, 2),mar=c(5,5,3,3))
  boxplot(NA,NA,NA,NA,
          boxwex=0.1,
          names=c("<70","70-180",">180","Total"),
          main = "Deviation Intervals", 
          xlab = "Comparator Glucose Range [mg/dL (mmol/L)]", 
          ylab = "Deviation [mg/dL or %]", 
          xaxs = "i",yaxs = "i",
          yaxt = "n", 
          xlim = c(0.635, 4.365), ylim = ylims)
  # Additional ticklabels
  xticks = c("(<3.9)","(3.9-10.0)","(>10.0)")
  for (x in 1:3){
    mtext(xticks[x],side=1,line=2,at=x)
  }
  # y-axis
  axis(2, at = seq(ylims[1], ylims[2], by = 20), las = 2, col = "black")
  
  # Background shading
  rect( 0.4, -15, 3.5,   15, col = adjustcolor("green4", 0.3), border = adjustcolor("green4", 0.3))
  rect( 0.4,  15, 3.5,   40, col = adjustcolor("darkorange", 0.3), border = adjustcolor("darkorange", 0.3))
  rect( 0.4, -15, 3.5,  -40, col = adjustcolor("darkorange", 0.3), border = adjustcolor("darkorange", 0.3))
  rect( 0.4,  40, 3.5, 1000, col = adjustcolor("red2", 0.4), border = 0)
  rect( 0.4, -40, 3.5,-1000, col = adjustcolor("red2", 0.4), border = 0)
  rect( 3.5, -20, 4.6,   20, col = adjustcolor("green4", 0.3), border = 0)
  rect( 3.5,  20, 4.6, 1000, col = adjustcolor("red2", 0.4), border = 0)
  rect( 3.5, -20, 4.6,-1000, col = adjustcolor("red2", 0.4), border = 0)
  
  # Horizontal lines
  abline(h = seq(ylims[1]+20, ylims[2]-20, by = 20 ), col = "gray60")
  abline(h = 0, col = "black")
  # Vertical Lines
  abline(v = c(1.5, 2.5, 3.5), col = "gray40", lty = 1, lwd = 2)
  
  # Boxes
  Int_size1 = c("85%","70%","80%","87%")
  Int_size2 = c("(98%)","(99%)","(99%)")
  
  wdth = 0.5
  blw = 1.5
  
  for (r in 1:4){
    
    TI_L1 = RES[r,"DI1_Lower"]
    TI_U1 = RES[r,"DI1_Upper"]
    TI_L2 = RES[r,"DI2_Lower"]
    TI_U2 = RES[r,"DI2_Upper"]
    
    if (r<=3) {
      # Interval 2
      rect(r-wdth/2,TI_L2,r+wdth/2,TI_U2,col=adjustcolor("grey60",0.4),
           border=adjustcolor("grey40",0.8),lwd=3.5)
      # Text
      mtext(Int_size2[r],line=0,side=3,at=r+0.21,cex=0.9,col="grey50")
    }
    
    # Interval 1
    rect(r-wdth/2,TI_L1,r+wdth/2,TI_U1,col=adjustcolor("grey50",0.8),
         border="black",lwd=3.5)
    # Text
    if (r<4){
      mtext(Int_size1[r],line=0,side=3,at=r-0.21,cex=0.9,col="black")
    }
    else {
      mtext(Int_size1[r],line=0,side=3,at=r,cex=0.9,col="black")
    }
    
    
    # Median
    med = RES[r,"Median"]
    rect(r-wdth/6,med,r+wdth/6,med,col="black",
         lwd=3.5)
    
  }
  
  # Version
  mtext(version,side=1,line=3.5,at=0,cex=0.6)
  
  # Secondary y-axis (mmol/L)
  par(new=TRUE)
  plot(1, 20 ,col="white",axes = FALSE, xlab = "", ylab = "",ylim=c(-80/18,80/18),yaxs = "i")
  axis(side=4,at=seq(-10, 10, by = 1 ), las=2)
  
  #--------------------------------------
  # Sensor-to-sensor variability
  # -------------------------------------------------------------------------
  
  # Setup axis
  par(mar=c(5,3,3,5))
  boxplot(NA,NA,NA,NA,
          boxwex=0.1,
          names=c("<70","70-180",">180","Total"),
          main = "Sensor-to-Sensor Variability",
          xlab = "Comparator Glucose Range [mg/dL (mmol/L)]",
          # ylab = "Deviation [mg/dL or %]",
          xaxs = "i",yaxs = "i",
          yaxt = "n",
          xlim = c(0.635, 4.365), ylim = c(-80, 80))
  # Additional ticklabels
  for (x in 1:3){
    mtext(xticks[x],side=1,line=2,at=x)
  }
  # y-axis
  axis(2, at = seq(ylims[1], ylims[2], by = 20), las = 2, col = "black")
  
  # Horizontal lines
  abline(h = seq(ylims[1]+20, ylims[2]-20, by = 20 ), col = "gray60")
  abline(h = 0, col = "black")
  # Vertical Lines
  abline(v = c(1.5, 2.5, 3.5), col = "gray40", lty = 1, lwd = 2)
  
  # Background
  rect( 0, -200, 5, 200, col = adjustcolor("gray80", 0.6), border = 0)
  
  # Errorbars
  
  # Number of sensors
  n_s <- length(unique(df$SensorID))
  
  # Get medians of sensors in total range
  ds_med <- subset(df, Range==4)
  ds_med <- data.frame(summarise(group_by(ds_med, SensorID),med=median(Diff)))
  ds_med <- ds_med[order(ds_med$med,decreasing=TRUE),]
  
  if (n_s>s_max){
    idx <- floor(linspace(1,n_s,s_max))
    ds_med <- ds_med[idx,]
    n_s <- s_max
  }
  
  x = seq(0.6,3.6,by=1)
  xi = linspace(0,0.8,n_s)
  
  for (r in 1:4){
    tmp <- subset(df,Range == r)
    for (i in 1:n_s){
      dev <- subset(tmp,SensorID == ds_med[i,"SensorID"])
      
      if (nrow(dev)>=3){
        m <- median(dev$Diff)
        if (nrow(dev)>=10){
          err_L = quantile(dev$Diff,probs=(1-dat_int)/2)
          err_U = quantile(dev$Diff,probs=(1+dat_int)/2)
          
          # Errorbars
          lines(x[r]+xi[i],m,type="p",pch=16,col="black")
          arrows(x0=x[r]+xi[i],x1=x[r]+xi[i],y0=err_L,y1=err_U,
                 code=3,angle=90,col="black",length=0)
        }
        else {
          err_L = quantile(dev$Diff,probs=0)
          err_U = quantile(dev$Diff,probs=1)
          
          # Errorbars
          lines(x[r]+xi[i],m,type="p",pch=16,col="black")
          arrows(x0=x[r]+xi[i],x1=x[r]+xi[i],y0=err_L,y1=err_U,
                 code=3,angle=90,col="black",length=0.015)
        }
      }
    }
  }
  
  # Secondary y-axis
  par(new=TRUE)
  plot(1, 20 ,axes = FALSE, xlab = "", ylab = "",ylim=c(-80/18,80/18),yaxs = "i")
  axis(side=4,at=seq(-10, 10, by = 1 ), las=2)
  mtext("Deviation [mmol/L]",side=4,line=3)
}



# ---------------------------------------
# Bootstrapping
bootstrapping <- function(df,N_BS,seed,conf_level=0.95,version="")
{
  
  # Calculate Acceleration --------------------------------------------------
  Calc_Acc <- function(df,qtl)
  {
    # Function to calculate the acceleration for BCa using a jackknife estimate with respect to the sensors
    # DiCiccio TJ, Efron B. Bootstrap confidence intervals. Stat Sci. 1996;11(3):189-228
    # Implementation is based on R package "bootstrap" (function bcanon)
    
    # Get set of sensors
    sens = unique(df$SensorID)
    n = length(sens)
    
    # Array with jackknife estimate
    # rows for sensors
    u = matrix(0,n,4)
    
    # Loop over sensors to remove one sensor at a time and re-estimate TIs
    for (i in 1:n){
      tmp = subset(df, SensorID != sens[i])
      u[i,] = quantile(tmp$Diff,probs=qtl)
    }
    
    # Estimate acceleration
    uu = sweep(-u,2,colMeans(u),FUN="+")
    a = colSums(uu^3) / (6 * (colSums(uu^2)^(3/2)))
    
    return(a)
  }
  
  # BCa ---------------------------------------------------------------------
  BCa <- function (dat, theta_h, a, qtl) 
  {
    # Function to calculate bias-corrected and accelerated bootstrap quantiles according to 
    # DiCiccio TJ, Efron B. Bootstrap confidence intervals. Stat Sci. 1996;11(3):189-228
    # Implementation is based on R package "bootstrap" (function bcanon)
    # In case of failure of the BCa method, the percentile method is used
    
    # Check if BCa can be applied
    if ((min(dat) >= theta_h) | (max(dat) <= theta_h)){
      # Switch to percentile method
      res = quantile(dat,probs=qtl)
      warning("BCa method could not be applied, using percentile method instead")
      return(res)
    }

    # BCa method
    sims <- length(dat)
    z.inv <- length(dat[dat < theta_h])/sims
    z <- qnorm(z.inv)
    
    # lower one-sided confidence interval
    bca_q <- pnorm(z + (z + qnorm(qtl))/(1 - a * (z + qnorm(qtl))))
    res <- quantile(dat, probs=bca_q)
    
    return(res)
  }
  
  
  # Boostrapping
  RES <- data.frame(Range=c("<70","70-180",">180","Total"),
                    Median=numeric(4),
                    DI1_Upper=numeric(4),
                    DI1_Lower=numeric(4),
                    DI1_Range=numeric(4),
                    DI2_Upper=numeric(4),
                    DI2_Lower=numeric(4),
                    DI2_Range=numeric(4))
  
  
  if (!(is.na(seed))){
    set.seed(seed)   
  }
  
  if (N_BS<10000){
    warning("Bootstrapping with less than 10 000 samples is not recommended")
  }
  
  Int_size1 = c(0.85,0.7,0.8,0.87)
  Int_size2 = c(0.98,0.99,0.99,0.87)
  
  # Loop over Ranges and calculate quantiles of intervals
  # Rows: Ranges
  # Columns: L1, U1, L2, U2
  DI = matrix(nrow=4,ncol=4)
  for (r in 1:4){
    tmp = subset(df,Range == r)
    DI[r,] = quantile(tmp$Diff,probs=c((1-Int_size1[r])/2,(1+Int_size1[r])/2,(1-Int_size2[r])/2,(1+Int_size2[r])/2))
    # Median
    RES[r,"Median"] = median(tmp$Diff)
  }
  
  BSdat <- as_tibble(df)
  BSdat <- BSdat %>% nest(data = -SensorID)
  
  BS <- bootstraps(BSdat, times = N_BS)
  
  stime <- Sys.time()
  print(paste("Bootstrapping N =",N_BS,"..."))
  # Loop over Ranges
  for (r in 1:4){
    res <- map(BS$splits, ~as_tibble(.) %>%
                 unnest(cols = c(data)) %>%
                 filter(Range == r) %>%
                 summarise(TI_L1=quantile(Diff,probs=(1-Int_size1[r])/2),
                           TI_U1=quantile(Diff,probs=(1+Int_size1[r])/2),
                           TI_L2=quantile(Diff,probs=(1-Int_size2[r])/2),
                           TI_U2=quantile(Diff,probs=(1+Int_size2[r])/2))) %>%
      bind_rows(.id = 'boots')
    
    # CI of AR using BCa
    # Get acceleration
    a <- Calc_Acc(subset(df, Range == r),c((1-Int_size1[r])/2,(1+Int_size1[r])/2,(1-Int_size2[r])/2,(1+Int_size2[r])/2))
    
    RES[r,"DI1_Lower"] = BCa(res$TI_L1,DI[r,1],a[1],(1-conf_level)/2)
    RES[r,"DI1_Upper"] = BCa(res$TI_U1,DI[r,2],a[2],(1+conf_level)/2)
    RES[r,"DI2_Lower"] = BCa(res$TI_L2,DI[r,3],a[3],(1-conf_level)/2)
    RES[r,"DI2_Upper"] = BCa(res$TI_U2,DI[r,4],a[4],(1+conf_level)/2)
    RES[r,"DI1_Range"] = RES[r,"DI1_Upper"] - RES[r,"DI1_Lower"]
    RES[r,"DI2_Range"] = RES[r,"DI2_Upper"] - RES[r,"DI2_Lower"]
    
    # Print progress
    print(paste("Range",r,"DONE"))
  }
  
  # Print Timing
  print(paste("Processing Time",round(difftime(Sys.time(),stime,units="secs"),2),"seconds"))
  
  # Delete results for Total range
  RES[4,"DI2_Lower"] = "-"
  RES[4,"DI2_Upper"] = "-"
  RES[4,"DI2_Range"] = "-"
  
  # Sensor to sensor variability parameters
  sens <- unique(df$SensorID)
  n_s <- length(unique(df$SensorID))
  
  for (r in 1:4){
    tmp <- subset(df,Range == r)
    med = numeric(n_s)
    r90 = numeric(n_s)
    
    for (i in 1:n_s){
      dev <- subset(tmp,SensorID == sens[i])
      
      if (nrow(dev)>=3){
        med[i] <- median(dev$Diff)
        if (nrow(dev)>=10){
          r90[i] = quantile(dev$Diff,probs=0.95) - quantile(dev$Diff,probs=0.05)
        }
        else {
          r90[i] = NA
        }
      }
      else{
        med[i] = NA
        r90[i] = NA
      }
    }
    
    RES[r,"BSV_Min_Max"] = paste("[",round(min(med,na.rm=TRUE),2)," - ",round(max(med,na.rm=TRUE),2),"]",sep="")
    RES[r,"BSV_Range"] = max(med,na.rm=TRUE) - min(med,na.rm=TRUE)
    RES[r,"WSV_90%Range"] = median(r90,na.rm=TRUE)
  }
  
  # Info
  RES[1,"Info"] = version
  RES[2,"Info"] = paste("N_BS:",N_BS,"Seed:",seed)
  RES[3,"Info"] = paste("Conf_Level:",conf_level)
  
  return (RES)
}



#-----------------------------------------
# Main function
CG_DIVA <- function(df,save_path,filename="CG-DIVA",
                    N_BS=10000,seed=1,
                    ylims=c(-80,80),s_max=25,figsize=c(30,15),
                    save_fig=TRUE,save_res=TRUE,show_fig=TRUE)
{
  
  version = "CG-DIVA v1.0"
  setDT(df) # Convert from data frame to data table
  
  ## Check Inputs
  # Columns
  for (col in c("SensorID","Comp","CGM")){
    if (!(col %in% colnames(df))){
      stop(paste("Column ",col," does not exist",sep=""))
    }
  }
  
  # NA
  if (any(is.na(df$SensorID)) | any(is.na(df$Comp)) | any(is.na(df$CGM))){
    stop("Dataset contains NA entries")
  }
  
  # Save Path
  if (!(dir.exists(save_path))){
    stop("Provided save_path does not exist")
  }
  
  ## Data processing 
  # calculate abs/rel deviations
  df$AbsDiff <- df$CGM - df$Comp
  df$RelDiff <- df$AbsDiff / df$Comp * 100
  # calculate TotDiff, i.e. abs/rel difference based on cutoff
  df <- df %>% mutate(Diff = if_else(Comp < 70, AbsDiff, RelDiff))
  # define concentration ranges
  df$Range <- if_else(df$Comp < 70, 1, 
                      if_else(df$Comp >= 70 & df$Comp <= 180, 2,
                              if_else(df$Comp > 180, 3, -1)))
  
  # Create new data table for total concentration range
  dfn <- df
  dfn$Range <- 4
  dfn$Diff <- df$RelDiff
  
  # Concatenate
  df <- rbind(df,dfn)
  
  for (r in 1:4){
    tmp = subset(df,Range == r)
    if (nrow(tmp)<100){
      stop(paste("Range",r,"contains an insufficient number of datapoints (<100)"))
    }
  }
    
  # bootstrapping
  RES <- bootstrapping(df,N_BS=N_BS,seed=seed,version=version)
  if (save_res){
    write.csv(RES,paste(save_path,filename,".csv",sep=""),row.names=FALSE)
  }
  

  # plotting
  if (save_fig){
    png(paste(save_path,filename,".png",sep=""),units="in",
        width=figsize[1]/2.5,height=figsize[2]/2.5,res=600)  
    plotting(df,RES,version,ylims=ylims,s_max=s_max,figsize=figsize)
    dev.off()
  }

  if (show_fig) {
    dev.new()
    plotting(df,RES,version,ylims=ylims,s_max=s_max,figsize=figsize)
  }  
  
}


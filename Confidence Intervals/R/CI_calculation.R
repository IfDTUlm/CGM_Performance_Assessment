# Calculation of confidence intervals on agreement rates for the assessment of compliance with FDA iCGM criteria

# For documentation see
# https://github.com/IfDTUlm/CGM_Performance_Assessment

# Created by: Institut für Diabetes-Technology Forschungs- und Entwichlungsgesellschaft mbH an der Universität Ulm
# Contact: m.eichenlaub@idt-ulm.de

# This is a free software and comes with ABSOLUTELY NO WARRANTY


# Load packages 
library(data.table)       # setDT and mutate
library(dplyr)            # %>% functionality
library(tidyr)            # nest functionality
library(purrr)            # map functionality
library(rsample)          # bootstrapping
library(GenBinomApps)     # Clopper-Pearson CI


# Clopper-Pearson Interval ------------------------------------------------
CI_Clopper_Pearson <- function (RES, df, alpha=0.05)
{
  # Calculate the Clopper-Pearson interval
  
  print("Calculate Clopper-Peason Intervals ...")

  n_k <- summarise(group_by(df, Range),
                   n=n(), 
                   k15=sum(WI15), 
                   k20=sum(WI20),
                   k40=sum(WI40))
  n_k <- data.frame(n_k)

  # Loop over Ranges
  for (r in 1:4){
    RES[r,"CP_CI15"] = clopper.pearson.ci(n_k[r,"k15"],n_k[r,"n"],alpha=alpha,CI="lower")$Lower.limit*100
    RES[r,"CP_CI20"] = clopper.pearson.ci(n_k[r,"k20"],n_k[r,"n"],alpha=alpha,CI="lower")$Lower.limit*100
    RES[r,"CP_CI40"] = clopper.pearson.ci(n_k[r,"k40"],n_k[r,"n"],alpha=alpha,CI="lower")$Lower.limit*100
  }
  
  return(RES)
}


# Wilson Interval ---------------------------------------------------------
CI_WilsonCC <- function (RES, df, alpha=0.05)
{
  # Calculate continuity corrected Wilson interval according to 
  # Short et al. "A novel confidence interval for a single proportion in the presence of clustered binary outcome data"
  # Statistical Methods in Medical Research 2020, Vol. 29(1) 111-121
  
  print("Calculate Wilson Intervals ...")

  lims = c(15,20,40)
  # Loop over Ranges
  for (r in 1:4){
    datr = subset(df, Range==r)
    # Loop over Limits
    for (i in 1:3){
      datr$WI = if_else(abs(datr$Diff) <= lims[i], 1, 0)
      
      datr = group_by(datr, SensorID)
      
      statr = summarize(datr,
                        m=n(),
                        x=sum(WI))
      m = statr$m
      x = statr$x
      
      n = nrow(statr)
      M1 = sum(m)
      M2 = sum(m^2)

      # Estimated AR
      pih = sum(x) / M1
      
      # Between cluster mean square
      BMS = (sum(x^2/m) - sum(x)^2 / M1) / (n - 1)
      
      # Within cluster mean square
      WMS = (sum(x) - sum(x^2/m)) / sum(m-1)
      
      # n*
      n_star = (M1^2 - M2) / ((n-1)*M1)
      
      # ICC
      rhoh = if_else(pih %in% c(0, 1), 1,
                     if_else(BMS - WMS < 0, 0,
                             (BMS - WMS) / (BMS + (n_star - 1) * WMS)))
      if (BMS - WMS < 0) {warning("ICC set to 0 due to negative correlation during Wilson interval calculation")}
      
      # VIF
      xih = 1 + rhoh*(M2-M1)/M1
      
      # WCC (lower one-sided confidence interval)
      z_a = qnorm(1 - alpha)
      WCC_CI = (2 * M1 * pih + xih * z_a^2 - 1 - sqrt(xih) * z_a * sqrt(xih * z_a^2 - 2 - 1 / M1 + 4 * pih * (M1 * (1 - pih) + 1))) / 
        (2 * (M1 + xih * z_a^2))
      
      # Collect Results
      if (i==1){ RES[r,"WCC_CI15"] = WCC_CI*100
      RES[r,"WCC_ICC15"] = rhoh}
      if (i==2){ RES[r,"WCC_CI20"] = WCC_CI*100
      RES[r,"WCC_ICC20"] = rhoh}
      if (i==3){ RES[r,"WCC_CI40"] = WCC_CI*100
      RES[r,"WCC_ICC40"] = rhoh}
      
    }
  }
  
  return(RES)
}


# Bootstrapping -----------------------------------------------------------
CI_Bootstrapping <- function (RES, df, alpha=0.05, N_BS=10000, seed=1)
{
  
  # Calculate Acceleration --------------------------------------------------
  Calc_Acc <- function(df)
  {
    # Function to calculate the acceleration for BCa using a jackknife estimate with respect to the sensors
    # DiCiccio TJ, Efron B. Bootstrap confidence intervals. Stat Sci. 1996;11(3):189-228
    # Implementation is based on R package "bootstrap" (function bcanon)
    
    # Get set of sensors
    subs = unique(df$SensorID)
    
    # Empty matrix, with columns for AR15, AR20 and AR40 and 
    # rows for sensors
    u = matrix(0,length(subs),3)
    
    # Loop over sensors to remove one sensor at a time and re-estimate ARs
    for (i in 1:length(subs)){
      tmp = subset(df, SensorID != subs[i])
      u[i,1] = sum(tmp$WI15)/length(tmp$WI15) * 100
      u[i,2] = sum(tmp$WI20)/length(tmp$WI20) * 100
      u[i,3] = sum(tmp$WI40)/length(tmp$WI40) * 100
    }
    
    # Estimate acceleration
    uu = sweep(-u,2,colMeans(u),FUN="+")
    a = colSums(uu^3) / (6 * (colSums(uu^2)^(3/2)))
    
    return(a)
  }
  
  # BCa ---------------------------------------------------------------------
  BCa <- function (dat, theta_h, a, alpha = 0.05) 
  {
    # Function to calculate bias-corrected and accelerated bootstrap quantiles according to 
    # DiCiccio TJ, Efron B. Bootstrap confidence intervals. Stat Sci. 1996;11(3):189-228
    # Implementation is based on R package "bootstrap" (function bcanon)
    
    sims <- length(dat)
    z.inv <- length(dat[dat < theta_h])/sims
    z <- qnorm(z.inv)
    
    # lower one-sided confidence interval
    prctl.inv <- pnorm(z + (z + qnorm(alpha))/(1 - a * (z + qnorm(alpha))))
    res <- quantile(dat, prctl.inv, names = FALSE, na.rm = TRUE)
    
    return(res)
  }
  
  if (!(is.na(seed))){
    set.seed(seed)   # for repeatability
  }
  
  if (N_BS<10000){
    warning("Bootstrapping with less than 10 000 samples is not recommended")
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
                 summarise(AR15=sum(WI15)/n()*100,
                           AR20=sum(WI20)/n()*100,
                           AR40=sum(WI40)/n()*100)) %>%
           bind_rows(.id = 'boots')
    
    # CI of AR using BCa
    # Get acceleration
    a <- Calc_Acc(subset(df, Range == r))
    
    RES[r,"BCa_CI15"] = BCa(res$AR15,RES[r,"AR15"],a[1],alpha=alpha)
    RES[r,"BCa_CI20"] = BCa(res$AR20,RES[r,"AR20"],a[2],alpha=alpha)
    RES[r,"BCa_CI40"] = BCa(res$AR40,RES[r,"AR40"],a[3],alpha=alpha)
    
    # Print progress
    print(paste("Range",r,"DONE"))
  }
  
  # Print Timing
  print(paste("Processing Time",round(difftime(Sys.time(),stime,units="secs"),2),"seconds"))
  
  # Collect Results
  RES[1,"Info"] = paste("Seed:",seed)
  RES[2,"Info"] = paste("N_BS:",N_BS)
  RES[3,"Info"] = paste("Conf_Level:",alpha)
  RES[4,"Info"] = paste("Software:",version[["version.string"]])
  
  
  return (RES)
}


# Main Function -----------------------------------------------------------
CI_Calculation <- function (df, save_path, filename="CI_Results",N_BS=10000, seed=1, alpha=0.05)
{
  # Calculation of confidence intervals on agreement rates using Clopper-Pearson, 
  # continuity-corrected Wilson and bootstrap methods
  
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
  df <- df %>% mutate(Diff = if_else(CGM < 70, AbsDiff, RelDiff))
  # define concentration ranges
  df$Range <- if_else(df$CGM < 70, 1, 
                          if_else(df$CGM >= 70 & df$CGM <= 180, 2,
                                  if_else(df$CGM > 180, 3, -1)))
  
  # Create new data table for total concentration range
  dfn <- df
  dfn$Range <- 4
  dfn$Diff <- df$RelDiff
  
  # Concatenate
  df <- rbind(df,dfn)
  
  # Calculate whether in or out of limits
  df <- df %>%
    mutate(WI15 = if_else(abs(Diff) <= 15, 1, 0),
           WI20 = if_else(abs(Diff) <= 20, 1, 0),
           WI40 = if_else(abs(Diff) <= 40, 1, 0))
  
  
  # Initialize results table
  RES <- data.frame(Range=c("<70","70-180",">180","Total"))
  
  # Calculate ARs
  res <- summarise(group_by(df, Range),
                   AR15=sum(WI15)/n()*100, 
                   AR20=sum(WI20)/n()*100, 
                   AR40=sum(WI40)/n()*100)
  
  RES[,"AR15"] = res$AR15
  RES[,"AR20"] = res$AR20
  RES[,"AR40"] = res$AR40
  
  
  ## Clopper-Pearson CI
  RES <- CI_Clopper_Pearson(RES,df,alpha=alpha)
  
  ## Wilson CC CI
  RES <- CI_WilsonCC(RES,df,alpha=alpha)
  
  ## Bootstrapping CI
  RES <- CI_Bootstrapping(RES, df, alpha=alpha, N_BS=N_BS, seed=seed)
  
  # Save Results 
  write.csv(RES,paste(save_path,filename,".csv",sep=""),row.names=FALSE)
  
  return (RES)
}
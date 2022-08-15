"""
Calculation of confidence intervals on agreement rates for the assessment of compliance with FDA iCGM criteria

For documentation see
https://github.com/IfDTUlm/CGM_Performance_Assessment

Created by: Institut für Diabetes-Technologie Forschungs- und Entwichlungsgesellschaft mbH an der Universität Ulm
Contact: cgm_performance@idt-ulm.de

This is a free software and comes with ABSOLUTELY NO WARRANTY
"""

import pandas as pd
import numpy as np
import scipy.stats as sts
import os
import warnings

def CI_Clopper_Pearson(RES, df, alpha=0.05):

    # Calculate the Clopper-Pearson interval
    
    print("Calculate Clopper-Pearson Intervals ...")

    n_k = df.groupby("Range")[["WI15","WI20","WI40"]].sum().reset_index()
    n_k["n"] = df.groupby("Range")["WI15"].count().to_numpy()

    # Loop over ranges
    for r in range(4):
        RES.at[r,"CP_CI15"] = sts.beta.ppf(alpha,n_k.at[r,"WI15"],
                                    n_k.at[r,"n"]-n_k.at[r,"WI15"]+1)*100
        RES.at[r,"CP_CI20"] = sts.beta.ppf(alpha,n_k.at[r,"WI20"],
                                    n_k.at[r,"n"]-n_k.at[r,"WI20"]+1)*100
        RES.at[r,"CP_CI40"] = sts.beta.ppf(alpha,n_k.at[r,"WI40"],
                                    n_k.at[r,"n"]-n_k.at[r,"WI40"]+1)*100

    return RES


def CI_WilsonCC(RES, df, alpha=0.05):

    # Calculate continuity corrected Wilson interval according to 
    # Short et al. "A novel confidence interval for a single proportion in the presence of clustered binary outcome data"
    # Statistical Methods in Medical Research 2020, Vol. 29(1) 111-121

    print("Calculate Wilson Intervals ...")
  
    # Loop over Ranges
    for r in range(4):
        datr = df[df["Range"] == r+1]
        # Loop over Limits
        for WI in ["WI15","WI20","WI40"]:
            
            statr = datr.groupby("SensorID")[WI].agg(["sum","count"]).reset_index()
            
            m = statr["count"].to_numpy()
            x = statr["sum"].to_numpy()            
      
            n = statr["SensorID"].nunique()
            M1 = np.sum(m)
            M2 = np.sum(m**2)
            
            # Estimated AR
            pih = np.sum(x) / M1
            
            # Between cluster mean square
            BMS = (np.sum(x**2/m) - np.sum(x)**2 / M1) / (n - 1)
            
            # Within cluster mean square
            WMS = (np.sum(x) - np.sum(x**2/m)) / np.sum(m-1)
            
            # n*
            n_star = (M1**2 - M2) / ((n-1)*M1)
            
            # ICC
            if pih == 0 or pih == 1:
                rhoh = 1
            else:
                if BMS-WMS < 0:
                    rhoh = 0
                    warnings.warn("ICC set to 0 due to negative correlation during Wilson interval calculation")
                else:
                    rhoh = (BMS - WMS) / (BMS + (n_star - 1) * WMS)           
            
            # VIF
            xih = 1 + rhoh * (M2 - M1) / M1
            
            # WCC (lower one-sided confidence interval)
            z_a = sts.norm.ppf(1-alpha)
            WCC_CI = (2 * M1 * pih + xih * z_a**2 - 1 - np.sqrt(xih) * z_a * np.sqrt(xih * z_a**2 - 2 - 1 / M1 + 4 * pih * (M1 * (1 - pih) + 1))) / (2 * (M1 + xih * z_a**2))
            
            # Collect Results
            RES.at[r,"WCC_CI"+WI[2:4]] = WCC_CI*100
            RES.at[r,"WCC_ICC"+WI[2:4]] = rhoh
            
    return RES


def CI_Bootstrapping(RES, df, alpha=0.05, N_BS=10000, seed=1):

    def calc_acc(df):
        """
        Function to calculate the acceleration for BCa using a jackknife estimate with respect to the sensors
        DiCiccio TJ, Efron B. Bootstrap confidence intervals. Stat Sci. 1996;11(3):189-228
        Implementation is based on R package "bootstrap" (function bcanon)

        Input:
        df:         Dataframe with data, reduced to a certain glucose range
        qtl:        List of quantiles for TIs  

        Output:
        a:      Numpy of acceleration for TI for provided quantiles
        

        """
        
        # Extract list of sensors
        sens = df["SensorID"].unique()
        sens.sort()
        n = len(sens)

        # Array with jackknife estimate
        u = np.zeros((n,3))    
        for i in range(n):
            # Remove single sensor and re-estimate 
            df_tmp = df[df["SensorID"] != sens[i]]
            u[i,:] = (df_tmp[["WI15","WI20","WI40"]].mean() * 100).to_numpy()
        
        # Remove mean
        uu = np.add(np.mean(u,axis=0),-u)
        # Estimate acceleration
        a = np.sum(uu**3,axis=0) / (6*(np.sum(uu**2,axis=0)**(3/2)))
        
        return a    
        
    def BCa(dat,theta_h,a,alpha=0.05):
        """
        Function to calculate bias-corrected and accelerated bootstrap quantiles according to 
        DiCiccio TJ, Efron B. Bootstrap confidence intervals. Stat Sci. 1996;11(3):189-228
        Implementation is based on R package "bootstrap" (function bcanon)
        In case of failure of the BCa method, the percentile method is used

        Input:
        dat:        Numpy array of bootstrapped samples
        theta_h:    Estimator of original sample
        a:          Acceleration
        alpha:      Quantile of boostrap sample to be calculated (0 to 1)

        Output:
        res:    Calculated quantile
        
        """
        
        # BCa method
        sims = len(dat)
        z_inv = np.sum((dat < theta_h)*1)/sims
        z = sts.norm.ppf(z_inv)

        prctl_inv = sts.norm.cdf(z + (z + sts.norm.ppf(alpha))/(1 - a * (z + sts.norm.ppf(alpha))))
        
        res = np.quantile(dat,prctl_inv)

        return res

    import time, sys

    ## Bootstrapping
    if seed:      # Seed is provided, if not reset
        np.random.seed(seed)       # For reproducibility
    else:
        np.random.seed()
    
    # Extract list of sensors
    sens = df["SensorID"].unique()
    sens.sort()
    
    # Empty Matrices containing bootstrap samples
    # Rows: Ranges, Columns: Interval limits/AR15,AR20,AR40, Depth: bootstrap samples 
    BS_AR = np.zeros((4,3,N_BS))
    
    # Define empty DataFrame for ARs so that missing values result in NaN instead of missing rows
    # This DataFrame has to have the same shape as what is intended to be converted to numpy
    df_emptyAR = pd.DataFrame(0, index=[0,1,2,3], columns=["WI15","WI20","WI40"])
    df_emptyAR["Range"] = [1,2,3,4]
    df_emptyAR.set_index("Range",inplace=True)

    # Start time for timing
    start = time.time()
    print("Bootstrapping (N = "+str(N_BS)+") ... ")
    # Loop over number of bootstap sample
    for i in range(N_BS):
        # Get clustered-bootstrap sample
        df_bs = pd.DataFrame({'SensorID':np.random.choice(sens, size=len(sens),replace=True)}).merge(df,how='left')
        
        # Agreement Rates
        # DataFrames are added to account for the case that a Range has no data (entry is NaN)
        BS_AR[:,:,i] = (df_emptyAR.add(df_bs.groupby("Range")[["WI15","WI20","WI40"]].mean()*100)).to_numpy()

        # Print status
        percent = int(np.round((i+1)/N_BS*100))
        bar = '+' * percent + "-" * (100-percent)
        print(f"\r|{bar}| {percent:.1f}%",end="\r")
    
    print()

    ## Collect Results
    # AR lower confidence bound
    # Loop over bins
    for r in range(4):
        # ARs
        a = calc_acc(df[df["Range"]==r+1])
        RES.at[r,"BCa_CI15"] = BCa(BS_AR[r,0,:],RES.at[r,"AR15"],a[0],alpha=alpha)
        RES.at[r,"BCa_CI20"] = BCa(BS_AR[r,1,:],RES.at[r,"AR20"],a[1],alpha=alpha)
        RES.at[r,"BCa_CI40"] = BCa(BS_AR[r,2,:],RES.at[r,"AR40"],a[2],alpha=alpha)

    # Print Timing    
    print("Processing Time: "+str(np.round(time.time()-start,2))+" seconds")

    # Provide Info
    RES.at[0,"Info"] = "Seed: "+str(seed)
    RES.at[1,"Info"] = "N_BS: "+str(N_BS)
    RES.at[2,"Info"] = "Conf_Level: "+str(alpha)

    return RES


def CI_calculation(df, save_path, filename="CI_Results",
                    N_BS=10000, seed=1, alpha=0.05):
 
   
    ## Check inputs
    # Dataframe columns
    for col in ["SensorID","Comp","CGM"]:
        if not(col in df.columns):
            raise ValueError("Column "+col+" does not exist")

    # Data Format
    if (pd.isna(df["SensorID"])).any() or (pd.isna(df["Comp"])).any() or (pd.isna(df["CGM"])).any():
        raise ValueError("Dataset contains NA entries")

    # Save path
    if not(os.path.isdir(save_path)):
        raise ValueError("Provided save_path does not exist")

    
    ## Data Processing
    df["AbsDiff"] = df["CGM"] - df["Comp"]                    # Absolute Difference
    df["RelDiff"] = df["AbsDiff"] / df["Comp"] * 100          # Relative Difference

    # Assign Ranges according to CGM
    # Range 1: <70, 2: 70-180, 3: >180
    df["Range"] = (df["CGM"]<70)*1 + ((df["CGM"]>=70) & (df["CGM"]<=180))*2 + (df["CGM"]>180)*3
    # Diff: Absolute Difference for <70, Relative Difference for >=70
    df["Diff"] = (df["CGM"]<70)*df["AbsDiff"] + (df["CGM"]>=70)*df["RelDiff"]

    # Making a copy of all data and assigning Ranges 4 (Total) and Relative difference to Diff
    dfn = df.copy()
    # Range 4: Total
    dfn["Range"] = np.ones(df.shape[0])*4
    # Diff: Relative difference
    dfn["Diff"] = df["RelDiff"]
    df = pd.concat([df,dfn])

    # Calculate whether in or out of limits
    df["WI15"] = (df["Diff"].abs() <= 15)*1
    df["WI20"] = (df["Diff"].abs() <= 20)*1
    df["WI40"] = (df["Diff"].abs() <= 40)*1

    # Initialize results table
    RES = pd.DataFrame() 
    RES["Range"] = ["<70","70-180",">180","Total"]

    # Calculate ARs
    RES["AR15"] = (df.groupby("Range")["WI15"].mean()*100).to_numpy()
    RES["AR20"] = (df.groupby("Range")["WI20"].mean()*100).to_numpy()
    RES["AR40"] = (df.groupby("Range")["WI40"].mean()*100).to_numpy()

    ## Clopper-Person CI
    RES = CI_Clopper_Pearson(RES,df,alpha=alpha)

    ## WilsonCC CI
    RES = CI_WilsonCC(RES, df, alpha=alpha)

    ## Bootstrapping CI
    RES = CI_Bootstrapping(RES, df, alpha=alpha, N_BS=N_BS, seed=seed)

    ## Save results
    RES.to_csv(save_path+filename+".csv",index=None)
    

"""
Continuous Glucose Deviation Interval and Variability Analysis (CG-DIVA)

For documentaiton see
https://github.com/IfDTUlm/CGM_Performance_Assessment

Created by: Institut für Diabetes-Technology Forschungs- und Entwichlungsgesellschaft mbH an der Universität Ulm
Contact cgm_performance@idt-ulm.de

This is a free software and comes with ABSOLUTELY NO WARRANTY
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

def boostrapping(df,N,seed,conf_level=0.95):

    def BCa(dat,theta_h,a,qtl):
        """
        Function to calculate bias-corrected and accelerated bootstrap quantiles according to 
        DiCiccio TJ, Efron B. Bootstrap confidence intervals. Stat Sci. 1996;11(3):189-228
        Implementation is based on R package "bootstrap" (function bcanon)
        In case of failure of the BCa method, the percentile method is used

        Input:
        dat:        Numpy array of bootstrapped samples
        theta_h:    Estimator of original sample
        a:          Acceleration
        qtl:        Quantile of boostrap sample to be calculated (0 to 1)

        Output:
        res:    Calculated quantile
        
        """
        
        import scipy.stats as sts

        # Remove NaNs
        dat = dat[~np.isnan(dat)]

        # Check if BCa can be applied
        if (np.min(dat) >= theta_h) or (np.max(dat) <= theta_h):
            # Switch to percentile method
            res = np.quantile(dat,qtl)
            return res

        # BCa method
        sims = len(dat)
        z_inv = np.sum((dat < theta_h)*1)/sims
        z = sts.norm.ppf(z_inv)

        if np.isinf(z):
            return np.NaN

        bca_q = sts.norm.cdf(z + (z + sts.norm.ppf(qtl))/(1 - a * (z + sts.norm.ppf(qtl))))
        
        res = np.quantile(dat,bca_q)

        return res


    def calc_acc(df,qtl=[]):
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
        u = np.zeros((n,4))    
        for i in range(n):
            # Remove single sensor and re-estimate 
            df_tmp = df[df["SensorID"] != sens[i]]
            u[i,:] = (df_tmp["Diff"].quantile(qtl)).to_numpy()
        
        # Remove mean
        uu = np.add(np.mean(u,axis=0),-u)
        # Estimate acceleration
        a = np.sum(uu**3,axis=0) / (6*(np.sum(uu**2,axis=0)**(3/2)))
        
        return a

    
    ## Bootstrapping
    import time

    if N<10000:
        print("WARNING: Bootstrapping with less than 10 000 samples is not recommended")

    RES = pd.DataFrame(columns=["Range","TI_U2","TI_U1","Median","TI_L1","TI_L2","Info"])
    
    if seed:      # Seed is provided, if not reset
        np.random.seed(seed)       # For reproducibility
    else:
        np.random.seed()
    
    # Interval Sizes according to FDA
    Int_size1 = [0.85,0.7,0.8,0.87]
    Int_size2 = [0.98,0.99,0.99,0.87]

    # Loop over Ranges and calculate quantiles of intervals
    # Rows: Ranges
    # Columns: L1, U1, L2, U2
    DI = np.zeros((4,4))
    for r in range(4):
        DI[r,:] = (df[df["Range"]==r+1]["Diff"].quantile([(1-Int_size1[r])/2,(1+Int_size1[r])/2,(1-Int_size2[r])/2,(1+Int_size2[r])/2])).to_numpy()

    # Extract list of sensors
    sens = df["SensorID"].unique()
    sens.sort()
    
    # Empty Matrices containing bootstrap samples
    # Rows: Ranges, Columns: Interval limits/AR15,AR20,AR40, Depth: bootstrap samples 
    BS_TI = np.zeros((4,4,N))
    
    # Start time for timing
    start = time.time()
    print("Bootstrapping (N = "+str(N)+") ... ")
    # Loop over number of bootstap sample
    for i in range(N):
        # Get clustered-bootstrap sample
        df_bs = pd.DataFrame({'SensorID':np.random.choice(sens, size=len(sens),replace=True)}).merge(df,how='left')
        
        # Loop over Ranges and calculate quantiles of intervals
        for r in range(4):
            BS_TI[r,:,i] = (df_bs[df_bs["Range"]==r+1]["Diff"].quantile([(1-Int_size1[r])/2,(1+Int_size1[r])/2,(1-Int_size2[r])/2,(1+Int_size2[r])/2])).to_numpy()
            
        # Print status
        percent = int(np.round((i+1)/N*100))
        bar = '+' * percent + "-" * (100-percent)
        print(f"\r|{bar}| {percent:.1f}%",end="\r")
    
    print()
    ## Collect Results
    RES["Range"] = ["<70","70-180",">180","Total"]
    RES.at[1,"Info"] = "N_BS: "+str(N)+" Seed: "+str(seed)
    RES.at[2,"Info"] = "Conf_Level: "+str(conf_level)

    # Median
    RES["Median"] = np.round((df.groupby("Range")["Diff"].median()).to_numpy(),2)

    # TI confidence intervals, AR lower confidence bound
    # Loop over Ranges
    for r in range(4):
        # TIs
        a = calc_acc(df[df["Range"]==r+1],qtl=[(1-Int_size1[r])/2,(1+Int_size1[r])/2,(1-Int_size2[r])/2,(1+Int_size2[r])/2])        
        RES.at[r,"TI_L1"] = np.round(BCa(BS_TI[r,0,:],DI[r,0],a[0],(1-conf_level)/2),2)
        RES.at[r,"TI_U1"] = np.round(BCa(BS_TI[r,1,:],DI[r,1],a[1],(1+conf_level)/2),2)
        RES.at[r,"TI_L2"] = np.round(BCa(BS_TI[r,2,:],DI[r,2],a[2],(1-conf_level)/2),2)
        RES.at[r,"TI_U2"] = np.round(BCa(BS_TI[r,3,:],DI[r,3],a[3],(1+conf_level)/2),2)

    # Delete results for Total range
    RES.at[3,"TI_L2"] = "-"
    RES.at[3,"TI_U2"] = "-"    
    
    # Print Timing    
    print("DONE","Processing Time: "+str(np.round(time.time()-start,2))+" seconds")
    
    return RES


def plotting(df,RES,version,ylims,s_max,figsize):

    def DI_plot(ax,RES,ylims=[-80,80]):
        """
        Plotting of Deviation Intervals in CG-DIVA

        Input:
        ax:         Handle to figure axis
        df:         DataFrame of results
        ylims:      Limits of y-axis
        
        """

        from matplotlib.patches import Rectangle as rect

        # Define Colors
        colors = [(94/255,156/255,32/255),      # Green    
                (251/255,145/255,36/255),       # Yellow   
                (238/255,12/255,4/255)]         # Red    
        al = 0.4    # Transparency of background colors

        ## Setup axis
        ax.set(xlim=[0.5,4.5],xticks=range(1,5),xticklabels=["<70\n(<3.9)","70-180\n(3.9-10.0)",">180\n(>10.0)","Total"],
            ylim=ylims,yticks=range(ylims[0],ylims[1]+1,20),
            xlabel="Comparator Glucose Range [mg/dL (mmol/L)]",
            ylabel="Deviation [mg/dL or %]")
        ax.set_title("Deviation Intervals",pad=15)     
        ax.grid(axis="y")
        ax2 = ax.twinx()
        ax2.set(ylim=np.array(ylims)/18)        

        # Plot vertical lines between boxes
        col, lw = 'grey', 1
        ax.plot([0,10],[0,0],color=col,linewidth=lw)
        for x in range(1,4):
            ax.plot([x+0.5]*2,ylims,color=col,linewidth=lw)

        # Colored Background    
        x1, x2 = [0,3.5], [3.5,4.5]
        # Green
        ax.fill_between(x1,[-15]*2,[15]*2,color=colors[0],alpha=al,edgecolor='none')
        ax.fill_between(x2,[-20]*2,[20]*2,color=colors[0],alpha=al,edgecolor='none')

        # Yellow
        ax.fill_between(x1,[15]*2,[40]*2,color=colors[1],alpha=al,edgecolor='none')
        ax.fill_between(x1,[-15]*2,[-40]*2,color=colors[1],alpha=al,edgecolor='none')

        # Red
        ax.fill_between(x1,[40]*2,[200]*2,color=colors[2],alpha=al,edgecolor='none')
        ax.fill_between(x1,[-40]*2,[-200]*2,color=colors[2],alpha=al,edgecolor='none')
        ax.fill_between(x2,[20]*2,[200]*2,color=colors[2],alpha=al,edgecolor='none')
        ax.fill_between(x2,[-20]*2,[-200]*2,color=colors[2],alpha=al,edgecolor='none')

        # Boxes of Deviation intervals
        wdth = 0.5  # Width of boxes
        blw = 1.5   # Line width of edges

        Int_size1 = ["85%","70%","80%","87%"]
        Int_size2 = ["(98%)","(99%)","(99%)"]

        # Loop over Ranges
        for r in range(4):    
            pos = r+1   # Position on x-axis

            # Read Results from DataFrame
            TI_L1, TI_U1 = RES.at[r,"TI_L1"], RES.at[r,"TI_U1"]
            TI_L2, TI_U2 = RES.at[r,"TI_L2"], RES.at[r,"TI_U2"]

            if r<3:         # Only plot Deviation interval 2 for ranges 1-3
                # Interval 2
                # Box            
                ax.add_patch(rect((pos-wdth/2,TI_L2),wdth,TI_U2-TI_L2,
                    facecolor="grey",alpha=0.4,edgecolor='k',linewidth=blw,linestyle="-",zorder=3))
                # Text above axis
                ax.text(pos+0.21,ax.get_ylim()[1]+2,Int_size2[r],ha="center",color="dimgrey")

            # Interval 1 
            # Box       
            ax.add_patch(rect((pos-wdth/2,TI_L1),wdth,TI_U1-TI_L1,
                facecolor="grey",alpha=0.8,edgecolor='k',linewidth=blw,linestyle="-",zorder=3))

            # Median
            Med = RES.at[r,"Median"]
            ax.plot([pos-wdth/2,pos+wdth/2],[Med]*2,'k-',linewidth=blw,zorder=3)

            # Text above axis
            if r<3:
                ax.text(pos-0.21,ax.get_ylim()[1]+2,Int_size1[r],ha="center")
            else:
                ax.text(pos,ax.get_ylim()[1]+2,Int_size1[r],ha="center")


    def S2SV_plot(ax,df,dat_int=0.9,ylims=[-80,80],s_max=25):
        """
        Plotting of Sensor-to-Sensor Variability in CG-DIVA

        Input:
        ax:         Handle to figure axis
        df:         DataFrame of raw data
        dat_int:    Interval of deviations to plot for each sensor (0 to 1)
        ylims:      Limits of y-axis
        
        """

        ## Setup axis
        ax.set(xlim=[0.5,4.5],xticks=range(1,5),xticklabels=["<70\n(<3.9)","70-180\n(3.9-10.0)",">180\n(>10.0)","Total"],
            ylim=ylims,yticks=range(ylims[0],ylims[1]+1,20),facecolor="lightgrey",
            xlabel="Comparator Glucose Range [mg/dL (mmol/L)]")
        ax.set_title("Sensor-to-Sensor Variability",pad=15)
        ax.grid(axis="y")
        ax2 = ax.twinx()
        ax2.set(ylim=np.array(ylims)/18,ylabel="Deviation [mmol/L]")

        # Plot vertical lines between boxes
        col, lw = 'grey', 1
        ax.plot([0,10],[0,0],color=col,linewidth=lw)
        for x in range(1,4):
            ax.plot([x+0.5]*2,ylims,color=col,linewidth=lw)

        # Number of sensors
        n_s = df["SensorID"].nunique()

        # Get medians of sensors of total range
        ds_med = df[df["Range"]==4].groupby("SensorID")["Diff"].median()
        ds_med = ds_med.reset_index()
        # Sort sensors (descending) according to median
        ds_med = ds_med.sort_values("Diff",ascending=False).reset_index(drop=True)
        if n_s > s_max:
            idx = np.floor(np.linspace(0,n_s-1,s_max)).astype(int)
            ds_med = ds_med.iloc[list(idx),:]         

        # Plotting
        x = np.arange(0.6,4.6)                   # Bin Position, plotting between 0.6 and 1.4, 1.6 and 2.4, ...
        xi = np.linspace(0,0.8,len(ds_med))      # Sensor Position within bin

        # Loop over Ranges
        for r in range(4):
            # Loop over Sensors
            for i,s in enumerate(ds_med["SensorID"].to_list()):
                # Get deviations of sensor in Range
                dev = (df[(df["SensorID"]==s) & (df["Range"]==r+1)]["Diff"]).to_numpy()

                if np.any(dev):         # Check if dev is empty
                    m = np.median(dev)
                    if len(dev)>=10:
                        err = [[m-np.quantile(dev,(1-dat_int)/2)],[np.quantile(dev,(1+dat_int)/2)-m]]
                        # Plot errorbar      
                        ax.errorbar(x[r]+xi[i],m,yerr=err,
                            marker='o',markersize=2,color='k',linewidth=0.75)
                        
                    else:   # Switch to full range of data
                        err = [[m-np.quantile(dev,0)],[np.quantile(dev,1)-m]]
                        # Plot errorbar with caps      
                        ax.errorbar(x[r]+xi[i],m,yerr=err,
                            marker='o',markersize=2,color='k',linewidth=0.75,capsize=1)
    
    # Plotting
    plt.rc('font',size=8)
    plt.rcParams.update({'font.sans-serif':'Arial'})
    
    print("Creating Figure ...")

    # Create figure
    fig, ax_a = plt.subplots(ncols=2,nrows=1,figsize=(figsize[0]/2.5,figsize[1]/2.5),constrained_layout=True)
    
    DI_plot(ax_a[0],RES,ylims=ylims)
    S2SV_plot(ax_a[1],df,ylims=ylims,s_max=s_max)

    # Print CG-DIVA info on plot
    ax_a[0].text(-0.2,-0.23,"CG-DIVA "+version,fontsize=6,
            transform=ax_a[0].transAxes)
    

def CG_DIVA(df,save_path,filename="CG-DIVA",
                N_BS=10000,seed=1,
                ylims=[-80,80],s_max=25,figsize=[16.5,8.5],
                save_fig=True,save_res=True,show_fig=True):
    """
    Perform CG-DIVA

    Inputs:
    df:         Pandas dataframe with columns "SensorID", "Comp" and "CGM"
    save_path:  Path for saving figure and results tables
    filename:   Filename of figure and results tables
    N_BS:       Number of samples for bootstrapping
    seed:       Seed for random number generator, provide [] when random seed shall be used
    ylims:      Limits of y-axis in CG-DIVA plot
    s_max:      Maximun number of sensor to display in sensor-to-sensor variability plot
    figsize:    [Width,Height] of figure
    save_fig:   True/False whether to save the figure
    save_res:   True/False whether to save the results in csv file
    show_plot:  True/False whether to show the figure 

    
    """
    version = "v1.0"
    print("\n\nCG_DIVA "+version)
    
    ## Check inputs
    # Dataframe columns
    for col in ["SensorID","Comp","CGM"]:
        if not(col in df.columns):
            raise ValueError("Column "+col+" does not exist")

    # Data Format
    if (pd.isna(df["SensorID"])).any():
        raise ValueError("Column SensorID contains NA entries")
    if df["Comp"].dtype != "float64" and df["Comp"].dtype != "int64":
        raise ValueError("Column Comp contains non-number entries")
    if df["CGM"].dtype != "float64" and df["CGM"].dtype != "int64":
        raise ValueError("Column CGM contains non-number entries")

    # Save path
    if not(os.path.isdir(save_path)):
        raise ValueError("Provided save_path does not exist")

    
    ## Data Processing
    df["AbsDiff"] = df["CGM"] - df["Comp"]                    # Absolute Difference
    df["RelDiff"] = df["AbsDiff"] / df["Comp"] * 100          # Relative Difference

    # Assign Ranges according to BG
    # Range 1: <70, 2: 70-180, 3: >180
    df["Range"] = (df["Comp"]<70)*1 + ((df["Comp"]>=70) & (df["Comp"]<=180))*2 + (df["Comp"]>180)*3 + df["Comp"]*0
    # Diff: Absolute Difference for <70, Relative Difference for >=70
    df["Diff"] = (df["Comp"]<70)*df["AbsDiff"] + (df["Comp"]>=70)*df["RelDiff"]

    # Making a copy of all data and assigning Ranges 4 (Total) and Relative difference to Diff
    dfn = df.copy()
    # Range 4: Total
    dfn["Range"] = np.ones(df.shape[0])*4
    # Diff: Relative difference
    dfn["Diff"] = df["RelDiff"]
    df = pd.concat([df,dfn])

    # Check data
    for r in range(4):
        n_r = df[df["Range"] == r+1]["Range"].count()
        if n_r<100:
            raise ValueError("Range "+str(r+1)+" CGM contains contains inssufficient number of datapoints (<100)")

    # Print dataset information
    print("Number of Datapoints:",int(df["Diff"].count()/2))
    print("Number of Sensors:",df["SensorID"].nunique())

    # Bootstrapping
    RES = boostrapping(df,N_BS,seed)
    RES.at[0,"Info"] = "CG-DIVA "+ version
    if save_res:
        RES.to_csv(save_path+filename+".csv",index=None)
    
    # Plotting
    
    plotting(df,RES,version,ylims,s_max,figsize)
    # Save and show plot
    if save_fig:
        plt.savefig(save_path+filename+".png",dpi=600)
    
    if show_fig:
        plt.show()
    print("\n\n")
    
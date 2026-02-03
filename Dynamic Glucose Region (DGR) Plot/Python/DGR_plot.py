"""
Dynamic Glucose Region plot

For documentaiton see
https://github.com/IfDTUlm/CGM_Performance_Assessment

Created by: Institut für Diabetes-Technology Forschungs- und Entwichlungsgesellschaft mbH an der Universität Ulm
Contact cgm_performance@idt-ulm.de

This is a free software and comes with ABSOLUTELY NO WARRANTY
"""
import matplotlib.pyplot as plt
import numpy as np
import os

# Define the parameters of the DGR plot
BGLow, BGHigh = 70, 300                     # Lower and upper BG limits
AlertLowBG, AlertHighBG = 80, 200           # Lower and upper BG limts for alert regions
AlertLowROC, AlertHighROC = -1, 1.5         # Lower and upper ROC limts for alert regions
pred_h = 30                                 # Prediction horizon for alert regions
ROC_lim = [-5,5]                            # Limits of ROC axis
BG_lim = [0,400]                            # Limits of BG axis
req = 7.5                                   # Minimum requirement for percentages in critical regions
pad = 0.2                                   # Padding between regions and bar plot (in units mg/dL/min)
bspace = 3.5                                # width of bar plot (in units mg/dl/min)


def region_cnt(df):
    """
    Count the number of RoC-BG pairs in each DGR plot region

    Input:
    df:     Dataframe with columns "BG" and "ROC"

    Output:
    cnt:    Array with number of pairs in each DGR region
            [BG low, BG high, Alert low, Alert high, Neutral, Total]

    Version History:

    """

    cnt = [0]*6
    n = df["BG"].count()

    # Predicted BG 30 min into the future
    df["BGP"] = df["BG"] + df["ROC"]*pred_h

    # BG low
    cnt[0] = (((df["BG"] < BGLow))*1).sum()
    # BG high
    cnt[1] = (((df["BG"] > BGHigh))*1).sum()
    # Alert low
    cnt[2] = (((df["BG"] >= BGLow) & (df["BGP"] < AlertLowBG) & (df["ROC"] < AlertLowROC))*1).sum()
    # Alert high
    cnt[3] = (((df["BG"] <= BGHigh) & (df["BGP"] > AlertHighBG) & (df["ROC"] > AlertHighROC))*1).sum()
    # Stable
    cnt[4] = (((df["BG"] <= 180) & (df["BG"] >= 70) & (df["ROC"] >= -1) & (df["ROC"] <= 1))*1).sum()
    # Total
    cnt[5] = n

    return np.array(cnt)


def remove_data(df):
    """
    Remove data in ellipitic region according to procedure
    described in the paper

    Input:
    df:         Dataframe with RoC-BG pairs

    Output:
    df:         Dataframe with RoC-BG pairs after exclusions
                (entries set to NaN)
    D_excl:     D values of exclusion ellipse
    msg:        Message string containing exclusion results
    n_ex:       Number of excluded datapoints

    Version History:

    """

    df = df.copy()
    df = df.dropna(subset=["BG","ROC"])

    cnt = region_cnt(df.copy())
    # Required number of datapoints in each region to fulfill the requirements
    req_cnt = np.ceil(req/100*cnt[5])

    # Check if requrements are fulfilled
    if all(cnt >= req_cnt):
        msg = "Data fulfill requirements"
        return df, 0, msg, 0
    else:
        # Number of datapoints that need to be excluded
        n_ex = int(np.ceil(cnt[5] - np.min(cnt) / (req/100)))
        # Shape parameters of ellipse
        a, b, BGc = 1, 115, 185
        # Calcualte D
        df["D"] = np.sqrt((df["ROC"]/a)**2+((df["BG"]-BGc)/b)**2)

        # Sort df according to D and create index for sorted df
        df = df.sort_values(by="D")
        # Identify the value of D for the exclusion ellipse
        D_excl = df.iloc[n_ex-1,2]

        if D_excl > 1:  # Ellipse extends into the critical regions
            msg = "Requirements cannot be fulfilled"
            df = df.sort_index()
            return df, D_excl, msg, 0

        else:
            # Insert NaNs for exluded pairs
            df.iloc[:n_ex,:] = np.nan
            # restore original sorting (I can't remember why)
            df = df.sort_index()
            msg = "{:d} ({:.1f}%) RoC-BG pairs excluded".format(n_ex,n_ex/cnt[5]*100)
            return df, D_excl, msg, n_ex


def DGR_plot(df,save_path=None,filename="DGR_plot",figsize=[13,10],ax=None,
             save_fig=False,show_fig=True,show_mmol=True,remove_dat=True):
    """
    Create the Dynamic Glucose Region plot based on BG-RoC pairs

    Inputs:
    df:             Pandas dataframe with columns "BG" and "ROC"
    save_path:      Path for saving figure and results
    filename:       Filename of figure and results
    figsize:        [Width,Height] of figure in cm
    ax:             Handle to existing axis object
    save_fig:       True/False whether to save the figure
    show_plot:      True/False whether to show the figure
    show_mmol:      True/False whether to inlcude axes in mmol/L or mmol/L/min
    remove_dat:     True/False whether to remove data to fullfill the requirements

    """

    ## Check inputs
    # Dataframe columns
    for col in ["BG","ROC"]:
        if not(col in df.columns):
            raise ValueError("Column "+col+" does not exist")

    # Data Format
    if df["BG"].dtype != "float64" and df["BG"].dtype != "int64":
        raise ValueError("Column BG contains non-number entries")
    if df["ROC"].dtype != "float64" and df["ROC"].dtype != "int64":
        raise ValueError("Column ROC contains non-number entries")

    # Save path
    if save_fig:
        if save_path is None:
            raise ValueError("Please provide save_path")
        else:
            if not(os.path.isdir(save_path)):
                raise ValueError("Provided save_path does not exist")

    # Remove data
    if remove_dat:
        df, _, msg, _ = remove_data(df)
        print(msg)

    #################
    # Setup the figure
    if ax is None:
        fig, ax = plt.subplots(figsize=[figsize[0]/2.5,figsize[1]/2.5],tight_layout=True)

    ax.set(xlim=[ROC_lim[0],ROC_lim[1]+bspace+pad],ylim=BG_lim,ylabel="BG concentration [mg/dL]",
            xticks=range(ROC_lim[0],ROC_lim[1]+1,1))
    ax.text(0.4,-0.15,"BG RoC [mg/dL/min]",ha="center",transform=ax.transAxes)

    if show_mmol:
        ax2_y = ax.twinx()
        ax2_y.set(ylim=[BG_lim[0]/18,BG_lim[1]/18],ylabel="BG concentration [mmol/L]",
                  yticks=range(int(BG_lim[0]/18),int(BG_lim[1]/18)+1,3))
        ax2_x = ax.twiny()
        ax2_x.set(xlim=[ROC_lim[0]/18,(ROC_lim[1]+bspace+pad)/18],
                  xticks=np.arange(-0.2,0.25,0.1))
        ax.text(0.4,1.11,"BG RoC [mmol/L/min]",ha="center",transform=ax.transAxes)

    # Axis Grid
    for x in np.arange(ROC_lim[0]+1,ROC_lim[1],1):
        ax.plot([x,x],BG_lim,color="grey",alpha=0.5,linewidth=0.5,zorder=0)
    for y in np.arange(BG_lim[0]+50,BG_lim[1],50):
        ax.plot(ROC_lim,[y,y],color="grey",alpha=0.5,linewidth=0.5,zorder=0)

    # Boundaries
    zo = 4
    border_lw=0.8
    col = "dimgrey"

    # BG low
    ax.plot(ROC_lim,[BGLow]*2,color=col,zorder=zo,linewidth=border_lw)
    # BG high
    ax.plot(ROC_lim,[BGHigh]*2,color=col,zorder=zo,linewidth=border_lw)
    # Alert low
    if AlertLowBG-ROC_lim[0]*pred_h <= BGHigh:
        ax.plot([AlertLowROC,AlertLowROC,ROC_lim[0]],
                [BGLow,AlertLowBG-AlertLowROC*pred_h,AlertLowBG-ROC_lim[0]*pred_h],color=col,zorder=zo,linewidth=border_lw)
    else:
        ax.plot([AlertLowROC,AlertLowROC,-(BGHigh-AlertLowBG)/pred_h],
                [BGLow,AlertLowBG-AlertLowROC*pred_h,BGHigh],color=col,zorder=zo,linewidth=border_lw)
    # Alert high
    if AlertHighBG-ROC_lim[1]*pred_h >= BGLow:
        ax.plot([AlertHighROC,AlertHighROC,ROC_lim[1]],
                [BGHigh,AlertHighBG-AlertHighROC*pred_h,AlertHighBG-ROC_lim[1]*pred_h],color=col,zorder=zo,linewidth=border_lw)
    else:
        ax.plot([AlertHighROC,AlertHighROC,-(BGLow-AlertHighBG)/pred_h],
                [BGHigh,AlertHighBG-AlertHighROC*pred_h,BGLow],color=col,zorder=zo,linewidth=border_lw)
    # Stable
    ax.plot([-1,1,1,-1,-1],[70,70,180,180,70],color=col,zorder=zo,linewidth=border_lw)

    # Borders between DGR and bars
    ax.plot([ROC_lim[1]]*2,BG_lim,color=col,zorder=zo,linewidth=1)
    ax.plot([ROC_lim[1]+pad]*2,BG_lim,color=col,zorder=zo,linewidth=1)

    # Region colors
    # BG high
    ax.fill([ROC_lim[0],ROC_lim[1],ROC_lim[1],ROC_lim[0]],
            [BG_lim[1],BG_lim[1],BGHigh,BGHigh],color="darkorange",alpha=0.4,zorder=2,edgecolor="None")
    # BG low
    ax.fill([ROC_lim[0],ROC_lim[1],ROC_lim[1],ROC_lim[0]],
            [BG_lim[0],BG_lim[0],BGLow,BGLow],color="darkorange",alpha=0.4,zorder=2,edgecolor="None")

    # Alert Low
    if AlertLowBG-ROC_lim[0]*pred_h <= BGHigh:
        ax.fill([ROC_lim[0],AlertLowROC,AlertLowROC,ROC_lim[0]],
                [BGLow,BGLow,AlertLowBG-AlertLowROC*pred_h,AlertLowBG-ROC_lim[0]*pred_h],
                color="red",alpha=0.2,zorder=2,edgecolor="None")
    else:
        ax.fill([ROC_lim[0],AlertLowROC,AlertLowROC,-(BGHigh-AlertLowBG)/pred_h,ROC_lim[0]],
                [BGLow,BGLow,AlertLowBG-AlertLowROC*pred_h,BGHigh,BGHigh],
                color="red",alpha=0.2,zorder=2,edgecolor="None")

    # Alert High
    if AlertHighBG-ROC_lim[1]*pred_h >= BGLow:
        ax.fill([ROC_lim[1],AlertHighROC,AlertHighROC,ROC_lim[1]],
                [BGHigh,BGHigh,AlertHighBG-AlertHighROC*pred_h,AlertHighBG-ROC_lim[1]*pred_h],
                color="red",alpha=0.2,zorder=2,edgecolor="None")
    else:
        ax.fill([ROC_lim[1],AlertHighROC,AlertHighROC,-(BGLow-AlertHighBG)/pred_h,ROC_lim[1]],
                [BGHigh,BGHigh,AlertHighBG-AlertHighROC*pred_h,BGLow,BGLow],
                color="red",alpha=0.2,zorder=2,edgecolor="None")

    # Neutral
    if AlertHighBG-ROC_lim[1]*pred_h >= BGLow:
        ax.fill([(AlertLowBG-BGHigh)/pred_h,AlertHighROC,AlertHighROC,ROC_lim[1],ROC_lim[1],AlertLowROC,AlertLowROC,(AlertLowBG-BGHigh)/pred_h],
                [BGHigh,BGHigh,AlertHighBG-AlertHighROC*pred_h,AlertHighBG-ROC_lim[1]*pred_h,BGLow,BGLow,AlertLowBG-AlertLowROC*pred_h,BGHigh],
                color="tab:green",alpha=0.2,zorder=2,edgecolor="None")
    else:
        ax.fill([(AlertLowBG-BGHigh)/pred_h,AlertHighROC,AlertHighROC,-(BGLow-AlertHighBG)/pred_h,AlertLowROC,AlertLowROC,(AlertLowBG-BGHigh)/pred_h],
                [BGHigh,BGHigh,AlertHighBG-AlertHighROC*pred_h,BGLow,BGLow,AlertLowBG-AlertLowROC*pred_h,BGHigh],
                color="tab:green",alpha=0.2,zorder=2,edgecolor="None")

    # Stable
    ax.fill([-1,1,1,-1],[70,70,180,180],color="tab:green",alpha=0.4,zorder=2,edgecolor="None")

    # Region labels
    alpha = 1
    pady = 10
    ax.text(ROC_lim[0]+0.2,350,"BG\nhigh",va="center",ha="left",color="orangered",alpha=alpha)
    ax.text(ROC_lim[1]-0.2,35,"BG\nlow",va="center",ha="right",color="orangered",alpha=alpha)
    ax.text(ROC_lim[1]-0.2,180,"Alert\nhigh",va="top",ha="right",color="red",alpha=alpha)
    ax.text(ROC_lim[0]+0.2,BGLow+pady,"Alert\nlow",va="bottom",ha="left",color="red",alpha=alpha)
    ax.text(1.8,60,"Stable",va="top",ha="center",color="tab:green",alpha=alpha)

    ##################################
    # Plot data
    # Scatter plot
    tmp = df.copy()
    tmp.loc[tmp["ROC"] > 5,"ROC"] = np.nan
    ax.scatter(tmp["ROC"],tmp["BG"],s=4,color="tab:blue",zorder=3,alpha=0.4,edgecolor="None")

    # Bars on the right
    cnt = region_cnt(df.copy())
    hmax = 12   # Maximum percentage of bar "axis"

    share = cnt / cnt[5] * 100
    share = [share[0],share[2],share[3],share[1],share[4]]

    # Bars for critical regions
    for i,col,alpha,pos,lbl,colt in zip(range(4),["darkorange","red","red","darkorange"],
                [0.4,0.2,0.2,0.4],[35,105,265,335],["BG low","Alert low","Alert high","BG high"],["orangered","red","red","orangered"]):
        ax.barh(pos,share[i]/hmax*bspace,zorder=3,color=col,alpha=alpha,height=45,left=ROC_lim[1]+pad)
        ax.text(ROC_lim[1]+pad*2,pos+10,lbl,va="center",ha="left",color=colt)
        ax.text(ROC_lim[1]+pad*2,pos-10,"{:.1f}%".format(share[i]),va="center",ha="left",color="black")

    # Number of pairs
    ax.text(ROC_lim[1]+pad*2,BG_lim[1]-5,"n={:d}".format(cnt[5]),va="top",ha="left",color="k")
    # Stable region
    ax.text(ROC_lim[1]+pad*2,195,"Stable",va="center",ha="left",color="tab:green")
    ax.text(ROC_lim[1]+pad*2,175,"{:.1f}%".format(share[4]),va="center",ha="left",color="black")
    # MARoC
    ax.text(ROC_lim[1]+bspace,BG_lim[0]-5,"MARoC: {:.2f}".format(df["ROC"].abs().mean()),va="top",ha="right",color="k")

    # Requirement
    ax.plot([ROC_lim[1]+pad+bspace*(req/hmax)]*2,[0,400],"--",color="red",linewidth=0.8)
    ax.text(ROC_lim[1]+pad+bspace*(req/hmax)+pad,BG_lim[1]-5,"{:.1f}%".format(req),color="red",va="top",ha="left")

    if save_fig:
        plt.savefig(save_path+filename+".png",dpi=600)

    if show_fig:
        plt.show()

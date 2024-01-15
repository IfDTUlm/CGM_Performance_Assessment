"""
Clincal Trend Concurrence Analysis

For documentaiton see
https://github.com/IfDTUlm/CGM_Performance_Assessment

Created by: Institut für Diabetes-Technology Forschungs- und Entwichlungsgesellschaft mbH an der Universität Ulm
Contact cgm_performance@idt-ulm.de

This is a free software and comes with ABSOLUTELY NO WARRANTY
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics as sklmetr
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle as rect
import os


def get_params(ncats):
    """
    Define basic parameters of the CTCA display

    Input:
    ncats:              Number of CGM RoC categories

    Output:
    ticklbls_comp:      xtick labels for comparator axis
    ticklbls_cgm:       xtick labels for cgm axis
    colors:             Matrix of colors (2D list) corresponding to the cells
    []:                 List of basic colors for zones A to D
    alpha:              alpha of colors

    Created:                ME 11.Dec 2023
    Checked/Versioned:      SP 11.Jan 2024

    Version History:

    """

    alpha = 0.8

    colA = (121/255,169/255,81/255)
    colB = (240/255,235/255,73/255)
    colC = (237/255,172/255,9/255)
    colD = (238/255,43/255,43/255)

    ticklbls_comp = ["<-3",
                     "-3 to <-2",
                     "-2 to <-1",
                     "-1 to +1",
                     ">+1 to +2",
                     ">+2 to +3",
                     ">+3"]

    if ncats == 7:
        ticklbls_cgm = ["(<-3) "+r"$\downdownarrows$",
                        "(-3 to <-2) "+r"$\downarrow$",
                        "(-2 to <-1) "+r"$\searrow$",
                        "(-1 to +1) "+r"$\rightarrow$",
                        "(>+1 to +2) "+r"$\nearrow$",
                        "(>+2 to +3) "+r"$\uparrow$",
                        "(>+3) "+r"$\upuparrows$"]

        colors = [[colA,colA,colB,colC,colD,colD,colD],
                [colA,colA,colA,colC,colD,colD,colD],
                [colB,colA,colA,colB,colC,colD,colD],
                [colD,colC,colB,colA,colB,colC,colC],
                [colD,colD,colC,colB,colA,colA,colB],
                [colD,colD,colD,colC,colA,colA,colA],
                [colD,colD,colD,colC,colB,colA,colA]]

    if ncats == 5:
        ticklbls_cgm = ["(<-2) "+r"$\downarrow$",
                        "(-2 to <-1) "+r"$\searrow$",
                        "(-1 to +1) "+r"$\rightarrow$",
                        "(>+1 to +2) "+r"$\nearrow$",
                        "(>+2) "+r"$\uparrow$"]

        colors = [[colA,colA,colA,colC,colD,colD,colD],
                  [colB,colA,colA,colB,colC,colD,colD],
                  [colD,colC,colB,colA,colB,colC,colC],
                  [colD,colD,colC,colB,colA,colA,colB],
                  [colD,colD,colD,colC,colA,colA,colA]]

    return ticklbls_comp, ticklbls_cgm, colors, [colA,colB,colC,colD], alpha


def zone_count(CMcnt,ncats):
    """
    Count the number of pairs in each zone

    Input:
    CMcnt:              Confusion matrix containing counts
    ncats:              Number of CGM RoC categories

    Output:
    Acnt:               Number of pairs in zone A
    Bcnt:               Number of pairs in zone B
    Ccnt:               Number of pairs in zone C
    Dcnt:               Number of pairs in zone D
    n:                  Total number of pairs

    Created:                ME 11.Dec 2023
    Checked/Versioned:      SP 11.Jan 2024

    Version History:

    """
    # The code sums of the respective cells for every zone row by row

    if ncats == 7:
        Acnt = 0
        for i,cols in enumerate([[0,1],[0,1,2],[1,2],[3],[4,5],[4,5,6],[5,6]]):
            Acnt = Acnt + CMcnt[i,cols].sum()
        Bcnt = 0
        for i,cols in enumerate([[2],[],[0,3],[2,4],[3,6],[],[4]]):
            Bcnt = Bcnt + CMcnt[i,cols].sum()
        Ccnt = 0
        for i,cols in enumerate([[3],[3],[4],[1,5,6],[2],[3],[3]]):
            Ccnt = Ccnt + CMcnt[i,cols].sum()
        Dcnt = 0
        for i,cols in enumerate([[4,5,6],[4,5,6],[5,6],[0],[0,1],[0,1,2],[0,1,2]]):
            Dcnt = Dcnt + CMcnt[i,cols].sum()
        n = np.sum(CMcnt)

    if ncats == 5:
        Acnt = 0
        for i,cols in enumerate([[0,1,2],[1,2],[3],[4,5],[4,5,6]]):
            Acnt = Acnt + CMcnt[i,cols].sum()
        Bcnt = 0
        for i,cols in enumerate([[],[0,3],[2,4],[3,6],[]]):
            Bcnt = Bcnt + CMcnt[i,cols].sum()
        Ccnt = 0
        for i,cols in enumerate([[3],[4],[1,5,6],[2],[3]]):
            Ccnt = Ccnt + CMcnt[i,cols].sum()
        Dcnt = 0
        for i,cols in enumerate([[4,5,6],[5,6],[0],[0,1],[0,1,2]]):
            Dcnt = Dcnt + CMcnt[i,cols].sum()
        n = np.sum(CMcnt)

    return Acnt, Bcnt, Ccnt, Dcnt, n


def CTCA(df,save_path=None,filename="CTCA",figsize=[16,6],fig=None,
             save_fig=False,show_fig=True):
    """
    Create the Dynamic Glucose Region plot based on BG-RoC pairs

    Inputs:
    df:             Pandas dataframe with columns "BG" and "ROC"
    save_path:      Path for saving figure and results
    filename:       Filename of figure and results
    figsize:        [Width,Height] of figure in cm
    fig:            Handle to existing figure object
    save_fig:       True/False whether to save the figure
    show_plot:      True/False whether to show the figure

    """

    ## Check inputs
    # Dataframe columns
    for col in ["Comp_ROC","Comp_ROC_Cat","CGM_ROC_Cat"]:
        if not(col in df.columns):
            raise ValueError("Column "+col+" does not exist")

    # Data Format
    if df["Comp_ROC"].dtype != "float64" and df["BG"].dtype != "int64":
        raise ValueError("Column Comp_ROC contains non-number entries")
    if df["Comp_ROC_Cat"].dtype != "float64" and df["BG"].dtype != "int64":
        raise ValueError("Column Comp_ROC_Cat contains non-number entries")
    if df["CGM_ROC_Cat"].dtype != "float64" and df["ROC"].dtype != "int64":
        raise ValueError("Column CGM_ROC_Cat contains non-number entries")

    # Save path
    if save_fig:
        if save_path is None:
            raise ValueError("Please provide save_path")
        else:
            if not(os.path.isdir(save_path)):
                raise ValueError("Provided save_path does not exist")

    ncats = df["Comp_ROC_Cat"].nunique()

    # Get parametes
    ticklbls_comp, ticklbls_cgm, colors, colors2, al = get_params(ncats)

    # Calculate Concurrence Matrix
    maroc = df["Comp_ROC"].abs().mean()
    tmp = df[["Comp_ROC_Cat","CGM_ROC_Cat"]].dropna()
    n = tmp.shape[0]
    # Confusion matrices need to be transposed to agree with the POCT layout (columns comparator, rows CGM)
    # Confusion matrix normalized over columns ("true" categories)
    CM = sklmetr.confusion_matrix(tmp["Comp_ROC_Cat"],tmp["CGM_ROC_Cat"],normalize="true").transpose()*100
    # Confusion matrix not normalized (just counts)
    CMcnt = sklmetr.confusion_matrix(tmp["Comp_ROC_Cat"],tmp["CGM_ROC_Cat"],normalize=None).transpose()
    # Remove top and bottom row for five-arrow system as they contain only zeros
    if ncats == 5:
        CM = CM[1:6,:]
        CMcnt = CMcnt[1:6,:]

    #################
    # Setup the figure
    if fig is None:
        fig = plt.figure(figsize=[figsize[0]/2.5,figsize[1]/2.5],tight_layout=True)

    gs = GridSpec(1,5,figure=fig)
    ax_a = [fig.add_subplot(gs[0,:4]),fig.add_subplot(gs[0,4])]

    # Setup axis
    ax = ax_a[0]
    ax.set(ylim=[-ncats-0.5,-0.5],xlim=[0.5,7.5],yticks=range(-ncats,0),
                yticklabels=reversed(ticklbls_cgm),xticks=[],xticklabels=[],
                ylabel="CGM RoC [mg/dl/min]")
    ax2 = ax.secondary_xaxis("top")
    ax2.set(xlim=[0.5,7.5],xticks=range(1,8),
            xticklabels=ticklbls_comp,xlabel="Comparator RoC [mg/dl/min]")

    # Plot lines
    yp = np.arange(-1,-ncats-1,-1)
    xp = np.arange(1,8,1)
    lw, col = 0.8, "black"
    for x in xp:
        ax.plot([x-0.5]*2,ax.get_ylim(),"-",linewidth=lw,color=col)
    for y in yp:
        ax.plot(ax.get_xlim(),[y+0.5]*2,"-",linewidth=lw,color=col)

    # Total
    txt_off = 0.03*(ax.get_ylim()[1]-ax.get_ylim()[0])
    ax.text(xp[0]-1,ax.get_ylim()[0]-txt_off,"Total",va="top",ha="center")

    # Plot colors and numbers in cells
    fs2 = 8
    off = 0.07
    # Loop over rows
    for i in range(ncats):
        # Loop over columns
        for j in range(7):
            # Column Total
            if i==0:
                ax.text(xp[j],ax.get_ylim()[0]-txt_off,"{:d}".format(CMcnt[:,j].sum()),
                            va="top",ha="center")
            # Color
            ax.add_patch(rect((xp[j]-0.5,yp[i]-0.5),1,1,
                        facecolor=colors[i][j],alpha=al,edgecolor=None))
            # Percentage number
            f_al = 1
            if CM[i,j] != 0:
                ax.text(xp[j],yp[i]-off,"{:.1f}%".format(CM[i,j]),color="k",
                        va="center",ha="center",fontsize=fs2,alpha=f_al)
            else:
                ax.text(xp[j],yp[i]-off,"-",
                        va="center",ha="center",fontsize=fs2,alpha=f_al)

    # Count Regions
    Acnt, Bcnt, Ccnt, Dcnt, n = zone_count(CMcnt,ncats)

    # Bar display on the right
    ax = ax_a[1]
    ax.axis("off")
    ax.set(ylim=[0,100],xlim=[0,1])
    bottom = 0
    xp = 0.25
    for cnt,col in zip([Acnt,Bcnt,Ccnt,Dcnt],colors2):
        ax.bar(xp,cnt/n*100,bottom=bottom,width=0.4,color=col,alpha=al)
        bottom = bottom + cnt/n*100

    # Percentages
    off = 0.25
    ax.text(xp+off,Acnt/n*100/2,"A: {:.1f}%".format(Acnt/n*100),va="center")
    ax.text(xp+off,(Acnt+Bcnt/2)/n*100,"B: {:.1f}%".format(Bcnt/n*100),va="center")
    ax.text(xp+off,(Acnt+Bcnt+Ccnt)/n*100,"C: {:.1f}%".format(Ccnt/n*100),va="top")
    ax.text(xp+off,100,"D: {:.1f}%".format(Dcnt/n*100),va="bottom")

    # MARoC
    ax.text(xp+off,0,"MARoC\n{:.2f}".format(maroc),va="bottom")

    # Total
    txt_off = 0.03*(ax.get_ylim()[1]-ax.get_ylim()[0])
    ax.text(xp,ax.get_ylim()[0]-txt_off,"{:d}".format(n),va="top",ha="center")


    if save_fig:
        plt.savefig(save_path+filename+".png",dpi=600)

    if show_fig:
        plt.show()


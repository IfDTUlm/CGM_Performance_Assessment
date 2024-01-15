Clinical Trend Concurrence Analysis (CTCA) as described in the following article (see pdf in the folder)

Eichenlaub, M.; Pleus, S.; Freckmann, G.: *A proposal for the clinical characterization of CGM trend arrow accuracy*, accepted the Journal of Diabetes Science and Technology

Should you decide to use this software in a publication we would appreciate if the above reference would be cited.

---

## Python
### Installation
All functions are included in the script *DGR_plot.py*.

Required packages
* pandas
* numpy
* matplotlib
* sklearn

### Usage

The main function is *CTCA* and has the following function call:

```
CTCA(df,save_path=None,filename="CTCA",figsize=[16,6],fig=None,
    save_fig=False,show_fig=True):
```
**Parameters:**

**df**: Pandas dataframe with columns *Comp_ROC*, *Comp_ROC_Cat* and *CGM_ROC_Cat*. *Comp_ROC* must contain the comparator ROCs in mg/dl/min. *Comp_ROC_Cat* must contain the categorized comparator RoCs typically encoded as -3,-2,-1,0,1,2,3. *CGM_ROC_Cat* must contain the categirized CGM ROC or documente CGM arrow, typicall as -3,-2,-1,0,1,2,3 for a seven-arrow and -2,-1,0,1,2 for a five-arrow CGM system.

**save_path**: Path for saving figure

**filename** *(optional)*: Filename of the output file *(default: CG_DIVA)*

**figsize** *(optional)*: [Width,Height] of saved figure in centimeters *(default [16.5,8.5])*

**ax** *(optional)*: Handle to existing figure object

**save_fig** *(optional)*: True/False whether to save the figure in png file *(default: True)*

**show_fig** *(optional)*: True/False whether to display the figure *(default: True)*

**Returns:**

Figure with the CG-DIVA plots as png file

An example of how to use the function is provided in the file *Example.py*.

### Example Figure

![](/Clinical Trend Concurrence Analysis/Python/CTCA.png)

## R

code to come ...


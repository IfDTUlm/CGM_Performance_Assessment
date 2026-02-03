Creation of a Dynamic Glucose Region (DGR) plot as described in the following articles

Eichenlaub, M.; Pleus, S.; Rothenbühler, M.; [...] Freckmann, G.: *Comparator Data Characteristics and Testing Procedures for the Clinical Performance Evaluation of Continuous Glucose Monitoring Systems*, Diabetes Technology & Therapeutics, 2024 [Link](https://www.liebertpub.com/doi/10.1089/dia.2023.0465)

Pleus, S.; Eichenlaub, M.; [...] Freckmann, G.: *Clinical assessment and acceptance criteria for continuous glucose monitoring (CGM) system performance: A proposed guideline by the IFCC Working Group on CGM*, Clinica Chimica Acta, 2025 [Link](https://www.sciencedirect.com/science/article/pii/S0009898125006072)

Should you decide to use this software in a publication we would appreciate if the above references would be cited.

## Update 2026
The function now inlcudes the updated definitions of the DGRs as proposed by the IFCC WG-CGM guideline published in 2025

---

## Python
### Installation
All functions are included in the script *DGR_plot.py*.

Required packages
* pandas
* numpy
* matplotlib

### Usage

The main function is *DGR_plot* and has the following function call:

```
DGR_plot(df,save_path=None,filename="DGR_plot",figsize=[13,10],ax=None,
             save_fig=False,show_fig=True,show_mmol=True,remove_dat=True):
```
**Parameters:**

**df**: Pandas dataframe with columns *BG* and *ROC*, containing the RoC-BG pairs

**save_path** *(optional)*: Path for saving the figure

**filename** *(optional)*: Filename of the output files *(default: DGR_plot)*

**figsize** *(optional)*: [Width,Height] of saved figure in centimeters *(default [13,10])*

**ax** *(optional)*: Handle to existing axis object

**save_fig** *(optional)*: True/False whether to save the figure in png file *(default: True)*

**show_fig** *(optional)*: True/False whether to display the figure *(default: True)*

**show_mmol** *(optional)*: True/False whether to include axes in mmol/L or mmol/L/min *(default: True)*

**remove_dat** *(optional)*: whether to remove data to fulfill the requirements *(default: True)*

**Returns:**

Figure with the CG-DIVA plots as png file and a csv file containing the deviation intervals

An example of how to use the function is provided in the file *Example.py*.

### Example Figure

![](</Dynamic Glucose Region (DGR) Plot/Python/DGR_plot.png>)

## R

available upon request


# CG-DIVA v1.0
Continuous Glucose Deviation Interval and Variability Analysis (CG-DIVA) as described in the following article

Eichenlaub, M.; Stephan, P.; Waldenmaier, D.; Pleus, S.; Rothenbühler, M.; Haug, C.; Hinzmann, R.; Thomas, A.; Jendle, J.; Diem, P.; Freckmann, G.: *Continuous Glucose Deviation Interval and Variability Analysis (CG-DIVA): A Novel Approach for the Statistical Accuracy Assessment of Continuous Glucose Monitoring Systems*, TBD 

Should you decide to use this software in a publication we would appreciate if the above reference would be cited.

---

## Python and R

### Installation

All functions are contained in one script that simply needs to be included.

Required packages for Python:

* pandas
* numpy
* matplotlib

Required packages for R:

* data.table
* dplyr
* tidyr
* purrr
* rsample
* pracma

### Usage

In both implementations, the main function *CG_DIVA* is called identically and performs the CG-DIVA. 

```
CG_DIVA(df,save_path,filename="CG-DIVA",
        N_BS=10000,seed=1,
        ylims=[-80,80],s_max=25,figsize=[16.5,8.5],
        save_fig=True,save_res=True,show_fig=True):
```
**Parameters:**

**df:** Pandas DataFrame (Python) or Data.Frame (R) with columns *SensorID*, *Comp* and *CGM*. *SensorID* must contain a unique identifier (string or number) for each CGM sensor. *Comp* must contain the  comparator measurement in mg/dL. *CGM* must contain the results of paired CGM measurements in mg/dL. Empty cells are not permitted

**save_path**: Path for saving figure and results tables

**filename** *(optional)*: Filename of the output files *(default: CG_DIVA)*

**N_BS** *(optional)*: Number of bootstrap samples for calculating the deviation intervals. For a reliable estimation of deviation intervals we recommend to use at least 10 000 bootstrap samples. However for testing purposes the sample size can be reduced to reduce computation time. *(default 10 000)*

**seed** *(optional)*: Seed for the random number generator used in the bootstrapping process. Provide [] if a the seed shall be automatically generated. Caution: Automatic seed generation can lead to slightly different results with each function call. To ensure reproducability provide a fixed seed *(default 1)*

**ylims** *(optional)*: Limits of y-axis (in mg/dL / %) in CG-DIVA figure *(default [-80,-80])*

**s_max** *(optional)*: Maximum number of sensors to be plotted in the sensor-to-sensor variability plot. If the number of sensors exeeds s_max, an equally spaced selection of sensors based on the median deviation in the total range (including maximum and minimum) are displayed *(default 25)*.

**figsize** *(optional)*: [Width,Height] of saved figure in centimeters

**save_fig** *(optional)*: True/False whether to save the figure in png file *(default: True)*

**save_res** *(optional)*: True/False whether to save the results of the deviation intervals in csv file *(default: True)*

**show_fig** *(optional)*: True/False whether to display the figure *(default: True)*

**Returns:**

Figure with the CG-DIVA plots as png file and a csv file containing the deviation intervals

An example of how to use the function is provided in the files *Example*. 

## Example Figures

### Python

![](/CG-DIVA/Python/CG-DIVA.png)

### R

![](/CG-DIVA/R/CG-DIVA.png)


# CG-DIVA v1.0
Continuous Glucose Deviation Interval and Variability Analysis (CG-DIVA) as described in the following article (see pdf in the folder)

XXX

Should you decide to use this software in a publication we would appreciate if the above reference would be cited.

---

## Python
### Installation
All functions are included in the script *CG_DIVA.py*. 

Required packages
* pandas
* numpy 
* matplotlib

### Usage

The main function is *CG_DIVA* and has the following function call:

```
CG_DIVA(df,save_path,filename="CG-DIVA",
        N_BS=10000,seed=1,
        ylims=[-80,80],s_max=25,figsize=[16.5,8.5],
        save_fig=True,save_res=True,show_fig=True):
```
**Parameters:**

**df**: Pandas dataframe with columns *SensorID*, *Comp* and *CGM*. *SensorID* must contain a unique identifier (string or number) for each CGM sensor. *Comp* must contain the  comparator measurement in mg/dL. *CGM* must contain the results of paired CGM measurements in mg/dL. Cells of *Comp* and *CGM* can be empty or contain *int* or *float* numbers. For the calculation of deviation intervals a minimum of 100 data pairs are required in each glucose range (<70, 70-180 and >180 mg/dL) 

**save_path**: Path for saving figure and results tables

**N_BS** *(optional)*: Number of bootstrap samples for calculating the deviation intervals. For a reliable estimation of deviation intervals we recommend to use at least 10 000 bootstrap samples. However for testing purposes the sample size can be reduced to reduce computation time. *(default 10 000)*

**seed** *(optional)*: Seed for the random number generator used in the bootstrapping process. Provide [] if a the seed shall be automatically generated. Caution: Automatic seed generation can lead to slightly different results with each function call. To ensure reproducability provide a fixed seed *(default 1)*

**ylims** *(optional)*: Limits of y-axis (in mg/dL / %) in CG-DIVA figure *(default [-80,-80])*

**s_max** *(optional)*: Maximum number of sensors to be plotted in the sensor-to-sensor variability plot. If the number of sensors exeeds s_max, an equally spaced selection of sensors based on the median deviation in the total range (including maximun and minimum) are displayed.

**figsize** *(optional)*: [Width,Height] of saved figure in centimeters *(default [16.5,8.5])*

**save_fig** *(optional)*: True/False whether to save the figure in png file *(default: True)*

**save_res** *(optional)*: True/False whether to save the results of the deviation intervals in csv file *(default: True)*

**show_fig** *(optional)*: True/False whether to display the figure *(default: True)*

**Returns:**

Figure with the CG-DIVA plots as png file and a csv file containing the deviation intervals

An example of how to use the function is provided in the file *Example.py*. 

### Example Figure

![](/CG-DIVA/Python/CG-DIVA_Test_Data.png)

## R

code to come ...


# Calculation of Confidence Intervals on Agreement Rates

Calculation of confidence intervals on agreement rates for the assessment of compliance with FDA iCGM criteria as described in 

XXReferenceXX

Should you decide to use this software in a publication we would appreciate if the above reference would be cited.

---

## Python and R

In both implementations, the main function *CI_calculation* is called identically and calculates the lower, one-sided 95% confidence intervals on the agreement rates as stipulated by the FDA using the Clopper-Pearson, clustered continuity-corrected Wilson and bias-corrected and accelerated bootstrap methods. 

```
CI_calculation(df,save_path,filename="CI_results",
                N_BS=10000,seed=1):
```
**Parameters:**

**df:** Pandas DataFrame (Python) or Data.Frame (R) with columns *SensorID*, *Comp* and *CGM*. *SensorID* must contain a unique identifier (string or number) for each CGM sensor. *Comp* must contain the  comparator measurement in mg/dL. *CGM* must contain the results of paired CGM measurements in mg/dL. Empty cells are not permitted

**save_path**: Path for saving results table

**filename** *(optional)*: Filename of the output files *(default: CI_results)*

**N_BS** *(optional)*: Number of bootstrap samples for calculating the bootstrapped confidence intervals. For a reliable estimation of we recommend to use at least 10 000 bootstrap samples. However for testing purposes the sample size can be reduced to reduce computation time. *(default 10 000)*

**seed** *(optional)*: Seed for the random number generator used in the bootstrapping process. Provide [] (Python) or NA (R) if a the seed shall be automatically generated. Caution: Automatic seed generation can lead to slightly different results with each function call. To ensure reproducability provide a fixed seed *(default 1)*

**Returns**:

A csv table with agreement rates (+/- 15 mg/dl or % (AR15), +/- 20 % (AR20), +/- 40 mg/dl or % (AR40)) and their lower, one-sided 95% confidence intervals as calculated by the three approaches.
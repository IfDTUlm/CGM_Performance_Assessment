import sys
sys.path.append("")     # Insert Path of CI_calculation.py
import CI_calculation as CI_calc
import pandas as pd

# Read data
df = pd.read_csv("") # Insert Path to data file

# Define save path
save_path = "" 

# Call function
CI_calc.CI_calculation(df,save_path)


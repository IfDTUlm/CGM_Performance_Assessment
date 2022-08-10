import sys
sys.path.append("")     # Insert Path of CG_DIVA.py
import CG_DIVA as CG_DIVA
import pandas as pd

# Read data
df = pd.read_csv("") # Insert Path to data file

# Define save path
save_path = "" 

# Perform CG-DIVA
CG_DIVA.CG_DIVA(df,save_path)

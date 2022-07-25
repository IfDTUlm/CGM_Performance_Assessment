import CG_DIVA as CG_DIVA
import pandas as pd

# Read Test data
df = pd.read_csv("CG-DIVA/Test_Data.csv")

# Define save path
save_path = "" 

# Perform CG-DIVA
CG_DIVA.CG_DIVA(df,save_path,filename="CG-DIVA_Test_Data")

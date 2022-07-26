import CI_calculation as CI_calc
import pandas as pd

# Read Test Data
df = pd.read_csv("CG-DIVA/Test_Data.csv")
# Define path to save results file
save_path = ""

# Call function
CI_calc.CI_calculation(df, save_path)


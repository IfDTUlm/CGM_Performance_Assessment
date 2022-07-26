import CI_calculation as CI_calc
import pandas as pd

df = pd.read_csv("V:/410-GM/IFCC - DCB/2137-DC/04-Auswertung/02-Publication CI/Software/Test_Data.csv")
save_path = "V:/410-GM/IFCC - DCB/2137-DC/04-Auswertung/02-Publication CI/Software/Python/"

CI_calc.CI_calculation(df, save_path)


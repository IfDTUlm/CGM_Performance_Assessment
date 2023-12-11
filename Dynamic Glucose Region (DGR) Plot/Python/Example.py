import pandas as pd
import matplotlib.pyplot as plt
import DGR_plot as dgr_plt

# Read Test data
df = pd.read_csv(".../Test_Data.csv")

# Define save path
save_path = ""

# Create DGR plot
plt.rcParams.update({"font.sans-serif":"Calibri","font.size":10})
dgr_plt.DGR_plot(df)

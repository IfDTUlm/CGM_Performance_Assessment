import pandas as pd
import matplotlib.pyplot as plt
import CTCA as ctca

# Read Test data
df = pd.read_csv("")

# Define save path
save_path = ""

# Create DGR plot
plt.rcParams.update({"font.sans-serif":"Arial","font.size":8.5})
ctca.CTCA(df)

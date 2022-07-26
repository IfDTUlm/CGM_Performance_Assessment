# Load Script
# Please provide the correct paths
source(".../CI_calculation.R")

# Please provide the correct paths
df <- read.csv(".../Test_Data.csv")
save_path <- ""

# Call Function
CI_Calculation(df,save_path)


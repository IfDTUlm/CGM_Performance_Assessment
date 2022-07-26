# Load Script
source("V:/410-GM/IFCC - DCB/2137-DC/04-Auswertung/02-Publication CI/Software/R/CI_calculation.R")

df <- read.csv("V:/410-GM/IFCC - DCB/2137-DC/04-Auswertung/02-Publication CI/Software/Test_Data.csv")
save_path <- "V:/410-GM/IFCC - DCB/2137-DC/04-Auswertung/02-Publication CI/Software/R/"

# Call Function
RES <- CI_Calculation(df,save_path)


import pandas as pd
import glob
import os
import matplotlib.pyplot as plt

# --- SETTINGS ---
folder = "weather/canada_historical_climate_data_timmins_airport_aug"
month  = 8
start_hour = 9     # Start of launch window
end_hour   = 17    # End of launch window inclusive

# --- LOAD ALL CSVs ---
csv_files = sorted(glob.glob(os.path.join(folder, "*.csv")))
print("Files found:", csv_files)

temps = []   # store all temperatures from all years

for file in csv_files:
    df = pd.read_csv(file)

    # Convert datetime
    df['Date'] = pd.to_datetime(df['Date/Time (LST)'])
    
    # Filter by month
    df = df[df['Date'].dt.month == month]

    # Extract hour
    df['Hour'] = df['Date'].dt.hour
    
    # Filter by hour window
    df = df[(df['Hour'] >= start_hour) & (df['Hour'] <= end_hour)]

    # Collect temperatures
    temps.extend(df['Temp (°C)'].dropna().tolist())

# --- PLOT HISTOGRAM ---
plt.figure(figsize=(10,6))
plt.hist(temps, bins=30, edgecolor='black', alpha=0.75)

plt.title(f"Histogram of Hourly Temperatures ({start_hour}:00–{end_hour}:00, Month={month})")
plt.xlabel("Temperature (°C)")
plt.ylabel("Frequency")
plt.grid(True, alpha=0.3)

plt.show()



import numpy as np
from math import erf, sqrt

temps = np.array(temps)

# Envelope limits
L, U = 17.0, 31.0

# Core statistics
mu = temps.mean()
sigma = temps.std(ddof=1)  # sample standard deviation

# z-scores
zL = (L - mu) / sigma
zU = (U - mu) / sigma

# Normal CDF (no scipy needed)
def norm_cdf(z):
    return 0.5 * (1 + erf(z / sqrt(2)))

# Coverage assuming normality
p_normal = norm_cdf(zU) - norm_cdf(zL)

# Empirical coverage
p_empirical = np.mean((temps >= L) & (temps <= U))

print("=== σ-based analysis ===")
print(f"Mean (μ)           = {mu:.2f} °C")
print(f"Std dev (σ)        = {sigma:.2f} °C")
print(f"Lower z (17°C)     = {zL:.2f}")
print(f"Upper z (31°C)     = {zU:.2f}")
print(f"P(17–31) normal    = {p_normal*100:.2f}%")
print(f"P(17–31) empirical = {p_empirical*100:.2f}%")


percentiles = [1, 2.5, 5, 10, 50, 90, 95, 97.5, 99]
values = np.percentile(temps, percentiles)

print("\n=== Percentiles of observed data ===")
for p, v in zip(percentiles, values):
    print(f"{p:5.1f}th percentile : {v:6.2f} °C")

# Where your limits sit
pct_below_L = np.mean(temps < L) * 100
pct_above_U = np.mean(temps > U) * 100

print("\n=== Envelope positioning ===")
print(f"% below 17°C = {pct_below_L:.2f}%")
print(f"% above 31°C = {pct_above_U:.2f}%")
print(f"% inside 17–31°C = {100 - pct_below_L - pct_above_U:.2f}%")

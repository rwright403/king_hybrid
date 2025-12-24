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

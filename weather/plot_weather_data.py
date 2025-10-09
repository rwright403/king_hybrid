import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Path to your CSVs
folder = "weather/canada_historical_climate_data_timmins_airport"

# Get all CSV files in the folder
csv_files = sorted(glob.glob(os.path.join(folder, "*.csv")))
print("Files found:", csv_files)

# Dictionaries for storing daily stats
morning_avg = {}   # 7:00–12:00
afternoon_avg = {} # 12:00–17:00

for file in csv_files:
    df = pd.read_csv(file)
    
    # Convert datetime
    df['Date'] = pd.to_datetime(df['Date/Time (LST)'])
    
    # Filter to August only
    august = df[df['Date'].dt.month == 8].copy()
    
    # Extract year from filename correctly
    year = os.path.basename(file).split('_')[5][-4:]
    
    # Extract hour
    august['Hour'] = august['Date'].dt.hour
    
    # Morning: 7–12
    morning = august[(august['Hour'] >= 7) & (august['Hour'] < 12)]
    morning_daily = morning.groupby(morning['Date'].dt.day)['Temp (°C)'].agg(['min', 'mean', 'max'])
    morning_avg[year] = morning_daily
    
    # Afternoon: 12–17
    afternoon = august[(august['Hour'] >= 12) & (august['Hour'] < 17)]
    afternoon_daily = afternoon.groupby(afternoon['Date'].dt.day)['Temp (°C)'].agg(['min', 'mean', 'max'])
    afternoon_avg[year] = afternoon_daily

# --- Morning plot ---
plt.figure(figsize=(10, 6))
for year, data in morning_avg.items():
    plt.plot(data.index, data['min'], label=f"{year} Min", linestyle='--', color='blue')
    plt.plot(data.index, data['mean'], label=f"{year} Mean", linestyle='--', color='green')
    plt.plot(data.index, data['max'], label=f"{year} Max", linestyle='--', color='red')

plt.title('Morning Temperatures (7am–12pm)')
plt.xlabel('Day of August')
plt.ylabel('Temperature (°C)')
plt.grid(True)
plt.legend(fontsize=8)
plt.show(block=False)  # <-- don't block

# --- Afternoon plot ---
plt.figure(figsize=(10, 6))
for year, data in afternoon_avg.items():
    plt.plot(data.index, data['min'], label=f"{year} Min", linestyle='--', color='blue')
    plt.plot(data.index, data['mean'], label=f"{year} Mean", linestyle='-', color='green')
    plt.plot(data.index, data['max'], label=f"{year} Max", linestyle='--', color='red')

plt.title('Afternoon Temperatures (12pm–5pm)')
plt.xlabel('Day of August')
plt.ylabel('Temperature (°C)')
plt.grid(True)
plt.legend(fontsize=8)
plt.show(block=False)  # <-- don't block


# Keep script alive to see both figures
input("Press Enter to exit...")
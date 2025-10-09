import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

"""
print(os.listdir("weather/canada_historical_climate_data_timmins_airport"))
"""

#NOTE: Path to your CSVs (change this!)
folder = "weather/canada_historical_climate_data_timmins_airport"

# Get all CSV files in the folder
csv_files = sorted(glob.glob(os.path.join(folder, "*.csv")))


print("Full glob path:", os.path.join(folder, "*.csv"))
print("Files found:", csv_files)

# Create 3 subplots
fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

# Loop through each CSV / year
for file in csv_files:
    df = pd.read_csv(file)
    
    # Convert 'Date/Time' to datetime
    df['Date'] = pd.to_datetime(df['Date/Time'])
    
    # Filter to only August
    august = df[df['Date'].dt.month == 8].copy()
    
    # Sort by day
    august.sort_values('Date', inplace=True)
    
    # Extract year from filename
    year = os.path.basename(file).split('.')[0]
    
    # Plot daily Min Temp
    axes[0].plot(august['Date'].dt.day, august['Min Temp (°C)'], marker='o', label=year)
    # Plot daily Mean Temp
    axes[1].plot(august['Date'].dt.day, august['Mean Temp (°C)'], marker='o', label=year)
    # Plot daily Max Temp
    axes[2].plot(august['Date'].dt.day, august['Max Temp (°C)'], marker='o', label=year)

# Configure subplots
axes[0].set_title('Daily Min Temperatures in August (All Years)')
axes[0].set_ylabel('°C')
axes[0].legend()
axes[0].grid(True)

axes[1].set_title('Daily Mean Temperatures in August (All Years)')
axes[1].set_ylabel('°C')
axes[1].legend()
axes[1].grid(True)

axes[2].set_title('Daily Max Temperatures in August (All Years)')
axes[2].set_xlabel('Day of August')
axes[2].set_ylabel('°C')
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.show()
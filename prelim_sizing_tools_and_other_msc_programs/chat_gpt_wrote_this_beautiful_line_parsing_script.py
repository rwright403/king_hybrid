import numpy as np
import matplotlib.pyplot as plt
import os
import cv2
import matplotlib
matplotlib.use('Agg')  # Use the Agg backend


# Load the image
img_path = os.path.join('prelim_sizing_tools_and_other_msc_programs', 'old_stuff_i_might_want_to_ref_later', 'tomasz_data.png')
#print("\n\n\n", img_path, "\nprelim_sizing_tools_and_other_msc_programs/old_stuff_i_might_want_to_ref_later/tomasz_data.png" "\n\n\n")
image = cv2.imread(img_path)


# Check if the image was loaded correctly
if image is None:
    print("Error: Could not read the image.")
else:
    print("Image loaded successfully.")
    print("Image shape:", image.shape)  # This will print (height, width, channels)

image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

# Define color ranges (Green, Orange, Blue) in RGB
green_lower = np.array([0, 120, 0])
green_upper = np.array([150, 255, 150])

orange_lower = np.array([150, 90, 0])
orange_upper = np.array([255, 180, 100])

blue_lower = np.array([0, 0, 100])
blue_upper = np.array([120, 120, 255])

# Function to create color masks
def create_mask(image, lower, upper):
    mask = cv2.inRange(image, lower, upper)
    return mask

# Create masks for each line
green_mask = create_mask(image_rgb, green_lower, green_upper)
orange_mask = create_mask(image_rgb, orange_lower, orange_upper)
blue_mask = create_mask(image_rgb, blue_lower, blue_upper)

# Extract line data from masks
def extract_line_data(mask):
    coords = np.column_stack(np.where(mask > 0))
    return coords

green_coords = extract_line_data(green_mask)
orange_coords = extract_line_data(orange_mask)
blue_coords = extract_line_data(blue_mask)

# Map pixel coordinates to graph coordinates
def map_to_graph_coords(coords, img_shape, x_min, x_max, y_min, y_max):
    x_pixels, y_pixels = img_shape[1], img_shape[0]
    x_mapped = x_min + (coords[:, 1] / x_pixels) * (x_max - x_min)
    y_mapped = y_max - (coords[:, 0] / y_pixels) * (y_max - y_min)
    return x_mapped, y_mapped


# Axis limits of your graph
time_min, time_max = 0, 14  # Time axis
pressure_min, pressure_max = 0, 800  # Pressure axis

# Convert pixel data to graph coordinates
green_time, green_pressure = map_to_graph_coords(green_coords, image.shape, time_min, time_max, pressure_min, pressure_max)
orange_time, orange_pressure = map_to_graph_coords(orange_coords, image.shape, time_min, time_max, pressure_min, pressure_max)
blue_time, blue_pressure = map_to_graph_coords(blue_coords, image.shape, time_min, time_max, pressure_min, pressure_max)

# Downsample to 100 points
def downsample(x, y, num_points=100):
    idx = np.round(np.linspace(0, len(x) - 1, num_points)).astype(int)
    return x[idx], y[idx]

green_time_ds, green_pressure_ds = downsample(green_time, green_pressure)
orange_time_ds, orange_pressure_ds = downsample(orange_time, orange_pressure)
blue_time_ds, blue_pressure_ds = downsample(blue_time, blue_pressure)

#convert pressure from psi to pa
def convert_pressure_to_si(color_pressure):
    si = []
    for p in color_pressure:
        p *= 1e5
        si.append(float(p))

    return si

green_pressure_ds = convert_pressure_to_si(green_pressure_ds)
orange_pressure_ds = convert_pressure_to_si(orange_pressure_ds)
blue_pressure_ds = convert_pressure_to_si(blue_pressure_ds)

# Plotting for verification
plt.figure(figsize=(10, 6))
plt.plot(green_time_ds, green_pressure_ds, label='Fuel Tank Pressure (Green)')
plt.plot(orange_time_ds, orange_pressure_ds, label='Oxidizer Tank Pressure (Orange)')
plt.plot(blue_time_ds, blue_pressure_ds, label='Combustion Chamber Pressure (Blue)')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.title('Static Fire 2 Pressure Data')
plt.legend()
plt.show()

# Save the data points
import pandas as pd

# Combine time and pressure into a single array, then sort by time
def sort_data_by_time(time, pressure):
    data = np.column_stack((time, pressure))
    sorted_data = data[data[:, 0].argsort()]  # Sort by the first column (time)
    return sorted_data[:, 0], sorted_data[:, 1]

# Sort each dataset by time
green_time_ds, green_pressure_ds = sort_data_by_time(green_time_ds, green_pressure_ds)
orange_time_ds, orange_pressure_ds = sort_data_by_time(orange_time_ds, orange_pressure_ds)
blue_time_ds, blue_pressure_ds = sort_data_by_time(blue_time_ds, blue_pressure_ds)

# Create DataFrames for each set of data
green_df = pd.DataFrame({'Time (s)': green_time_ds, 'Pressure (psi)': green_pressure_ds})
orange_df = pd.DataFrame({'Time (s)': orange_time_ds, 'Pressure (psi)': orange_pressure_ds})
blue_df = pd.DataFrame({'Time (s)': blue_time_ds, 'Pressure (psi)': blue_pressure_ds})

# Save each DataFrame to a CSV file
green_df.to_csv('green_fuel_tank_pressure.csv', index=False)
orange_df.to_csv('orange_oxidizer_tank_pressure.csv', index=False)
blue_df.to_csv('blue_combustion_chamber_pressure.csv', index=False)

print("Data saved to CSV files!")
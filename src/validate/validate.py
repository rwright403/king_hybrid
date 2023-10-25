import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Read the CSV file
file_path = 'your_data.csv'  # Replace with the path to your CSV file
modeldf = pd.read_csv('./src/thrust.csv', header=None)

# Step 2: Prepare data
# Assuming your CSV has columns 'x' and 'y', change these column names accordingly
x = modeldf.iloc[:, 0]
y = modeldf.iloc[:, 1]

# Step 3: Create a simple line graph
plt.figure(figsize=(10, 6))  # Optional: Adjust the figure size
plt.plot(x, y, marker='o', linestyle='-')

# Step 4: Customize the graph (optional)
plt.title('CSV Data Visualization')
plt.xlabel('X-axis Label')
plt.ylabel('Y-axis Label')

# Step 5: Show the graph or save it to a file
# To display the graph:
plt.show()

# To save the graph as an image (e.g., PNG):
# plt.savefig('output_graph.png')

# Optional: Close the plot window if you're not displaying the graph interactively
# plt.close()

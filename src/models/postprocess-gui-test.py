"""
import pyvista as pv
import numpy as np

# Create a dummy FEA mesh (replace with your actual node coordinates & connectivity)
nodes = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1]])
cells = [4, 0, 1, 2, 3]  # Quad element (modify for your element type)

# Create the mesh
grid = pv.PolyData(nodes)
grid["stress"] = np.random.rand(len(nodes))  # Example stress values

# Plot with stress color map
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars="stress", point_size=10, cmap="coolwarm")
plotter.show()
"""

import pyvista as pv
import numpy as np
import time

# Create a simple FEA mesh (replace with your actual node coordinates)
nodes = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1]])

# Define the elements (can be tetrahedral, hexahedral, etc. based on your mesh)
cells = [4, 0, 1, 2, 3]  # A single tetrahedron element (modify as needed)

# Create PyVista mesh
grid = pv.PolyData(nodes)
grid.faces = np.array(cells)

# Dummy time-varying stress data (replace with actual stress values)
num_timesteps = 50
stress_over_time = np.linspace(0, 1, num_timesteps)[:, None] * np.random.rand(len(nodes))

# Dummy displacement data: apply time-varying factor to random displacements for each node
displacements = 0.1 * np.sin(np.linspace(0, 2 * np.pi, num_timesteps))[:, None]  # Shape (50, 1)
random_displacements = np.random.rand(len(nodes), 3)  # Random displacements for each node, shape (5, 3)
displacements = displacements[:, None, :] * random_displacements  # Broadcasting to apply displacement across nodes

# Setup PyVista plotter for animation
plotter = pv.Plotter(off_screen=False)  # Ensure off_screen is False for interaction
actor = plotter.add_mesh(grid, scalars=stress_over_time[0], point_size=10, cmap="coolwarm")

# Initialize the render window and interactor (this is necessary before rendering)
plotter.show(auto_close=False)

# Animate the stress and displacement change over time
for t in range(num_timesteps):
    # Update displacement for node motion
    grid.points = nodes + displacements[t]
    
    # Update the stress as a scalar array for color mapping
    grid.point_arrays["stress"] = stress_over_time[t]  # Assigning to point arrays
    
    # Update the mesh in the plotter
    plotter.update()

    # Render the scene (display)
    plotter.render()
    
    # Optional: Add a delay for smooth animation
    time.sleep(0.05)

# Show the final animation
plotter.show()

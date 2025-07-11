from paraview.simple import *
import os

# Set input/output directories
vtk_dir = "simu1D/output_lat_30/vtk_pngs/slice1"
output_dir = "simu1D/output_lat_30/vtk_pngs/slice1"
os.makedirs(output_dir, exist_ok=True)

# Get and sort .vtk files
vtk_files = sorted(f for f in os.listdir(vtk_dir) if f.endswith(".vtk"))

# Load first file to set up camera/view
reader = LegacyVTKReader(FileNames=[os.path.join(vtk_dir, vtk_files[0])])
display = Show(reader)
view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [800, 600]
view.Background = [1, 1, 1]  # white

# Optional: color by scalar (skip or modify as needed)
ColorBy(display, ('POINTS', 'your_scalar_name'))

# Set camera (optional: adjust to fit your data)
ResetCamera()
Render()

# Loop through and save images
for i, filename in enumerate(vtk_files):
    reader.FileNames = [os.path.join(vtk_dir, filename)]
    Render()
    SaveScreenshot(os.path.join(output_dir, f"frame_{i:04d}.png"), view)
    print(f"Saved frame_{i:04d}.png")


#
reader = LegacyVTKReader(FileNames=["simu1D/output_lat_30/vtk_pngs/slice1/wave100.vtk"])
Show(reader)
Render()
data = servermanager.Fetch(reader)

# List point and cell data arrays
print("Point data arrays:", [data.GetPointData().GetArrayName(i) for i in range(data.GetPointData().GetNumberOfArrays())])
print("Cell data arrays:", [data.GetCellData().GetArrayName(i) for i in range(data.GetCellData().GetNumberOfArrays())])

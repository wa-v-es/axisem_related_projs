# run this through terminal
# pvbatch render_vtk_png.py
from paraview.simple import *
import os
import sys

reader = LegacyVTKReader(FileNames=["output_10sec/elements/orthogonal_azimuthal_slices/vtk/slice1/wave100.vtk"])
Show(reader)
Render()
data = servermanager.Fetch(reader)
# print('len of data',len(data))
# print('len of data',len(data.GetPointData().GetNumberOfArrays()))

# List point and cell data arrays
print("Point data arrays:", [data.GetPointData().GetArrayName(i) for i in range(data.GetPointData().GetNumberOfArrays())])
print("Cell data arrays:", [data.GetCellData().GetArrayName(i) for i in range(data.GetCellData().GetNumberOfArrays())])

# sys.exit()
# Set input/output directories
vtk_dir = "output_10sec/elements/orthogonal_azimuthal_slices/vtk/slice1"
output_dir = "output_10sec/elements/orthogonal_azimuthal_slices/vtk_pngs/slice1"
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
ColorBy(display, ('POINTS', 'UZ'))

# Set camera (optional: adjust to fit your data)
ResetCamera()
Render()

# Loop through and save images
for i, filename in enumerate(vtk_files):
    reader.FileNames = [os.path.join(vtk_dir, filename)]
    Render()
    SaveScreenshot(os.path.join(output_dir, f"frame_{i:04d}.png"), view)
    print(f"Saved frame_{i:04d}.png")

sys.exit()
#

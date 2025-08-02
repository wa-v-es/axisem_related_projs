
#pvpython
from paraview.simple import *
import os
import sys

# Filename stuff
vtk_dir = "../output_10sec/elements/orthogonal_azimuthal_slices/vtk/slice3/"
output_dir = "../output_10sec/elements/orthogonal_azimuthal_slices/vtk_pngs/slice3/"
os.makedirs(output_dir, exist_ok=True)
# path_to_vtk  = '/Users/eaton/Downloads/'
wave_id      = 100
slicename = f'wave{wave_id}'

# Load the VTK file
vtk_files = sorted([f for f in os.listdir(vtk_dir) if f.endswith('.vtk')])
for i, vtk_file in enumerate(vtk_files):
# for vtk_file in vtk_files:
    outpath = os.path.join(output_dir, f"{vtk_file[:-3]}png")
    wave100vtk = LegacyVTKReader(FileNames=[os.path.join(vtk_dir, vtk_file)])
    # wave100vtk = LegacyVTKReader(registrationName=f'{slicename}.vtk',
    #                              FileNames=[f'{path_to_vtk}/{slicename}.vtk'])

    # Get the 'View' - this may cause an issue from the terminal??
    rv = GetActiveViewOrCreate('RenderView')
    SetActiveView(rv)

    # Display the vtk data
    wave100vtkDisplay = Show(wave100vtk, rv, 'UnstructuredGridRepresentation')

    # Adjust the camera to be orthogonal to slice
    # rv.ResetActiveCameraToPositiveY()# slice0
    rv.ResetActiveCameraToPositiveX()# slice0

    rv.ResetCamera(False, 0.9) # slice0
    wave100vtkDisplay = Show(wave100vtk, rv, 'UnstructuredGridRepresentation')

    # Change colourbar limits
    vmax = 1e-8
    vmin = -vmax
    uZLUT = GetColorTransferFunction('UZ')
    uZLUT.RescaleTransferFunction(vmin, vmax)

    # Render (update) before taking screenshot
    # Save with transparent background
    SaveScreenshot(filename=outpath,
                   viewOrLayout=rv,
                   ImageResolution=[1050, 543],#[2100, 1086]
                   TransparentBackground=1)
    print("Saved to ", outpath)
    # sys.exit()

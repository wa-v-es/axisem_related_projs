#magick -delay 80 -quiet -quality 50 -loop 0 *.png out.gif
#pvpython
from paraview.simple import *
import os
import sys
import re

# Filename stuff
vtk_dir = "../plumes_iaspi91_10sec_new_loc_with_wave/simu1D/output/elements/orthogonal_azimuthal_slices/vtk/slice3/"
output_dir = "../plumes_iaspi91_10sec_new_loc_with_wave/simu1D/output/elements/orthogonal_azimuthal_slices/vtk_pngs/slice3/"
os.makedirs(output_dir, exist_ok=True)
# path_to_vtk  = '/Users/eaton/Downloads/'
wave_id      = 100
slicename = f'wave{wave_id}'

# Load the VTK file
vtk_files = sorted([f for f in os.listdir(vtk_dir) if f.endswith('.vtk')])
print(vtk_files[::15])
# sys.exit()
for i, vtk_file in enumerate(vtk_files[::15]):
# for vtk_file in vtk_files:
    match = re.search(r"wave(\d+)\.vtk", vtk_file)
    if match:
        num = int(match.group(1))
        padded_name = f"wave{num:04d}.png"  # Pads to 4 digits: 0001, 0190, etc.
    else:
        padded_name = f"{vtk_file[:-3]}.png"

    outpath = os.path.join(output_dir, padded_name)
    # outpath = os.path.join(output_dir, f"{vtk_file[:-3]}png")
    wave100vtk = LegacyVTKReader(FileNames=[os.path.join(vtk_dir, vtk_file)])
    # wave100vtk = LegacyVTKReader(registrationName=f'{slicename}.vtk',
    #                              FileNames=[f'{path_to_vtk}/{slicename}.vtk'])

    # Get the 'View' - this may cause an issue from the terminal??
    rv = GetActiveViewOrCreate('RenderView')
    SetActiveView(rv)

    # Display the vtk data
    wave100vtkDisplay = Show(wave100vtk, rv, 'UnstructuredGridRepresentation')

    # Adjust the camera to be orthogonal to slice
    # rv.ResetActiveCameraToPositiveY()# slice0,2
    rv.ResetActiveCameraToPositiveX()# slice1,3

    rv.ResetCamera(False, 0.9) # slice0
    wave100vtkDisplay = Show(wave100vtk, rv, 'UnstructuredGridRepresentation')

    # Change colourbar limits
    vmax = 5e-9
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

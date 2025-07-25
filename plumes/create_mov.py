import os
import vtk
import cv2
import numpy as np
from vtk.util import numpy_support

# Directory containing the .vtk files
vtk_dir = "simu1D/output/elements/orthogonal_azimuthal_slices/vtk/slice2/"
output_video = "output_slice2.mov"
frame_rate = 5  # frames per second

# Get sorted list of .vtk files
vtk_files = sorted([f for f in os.listdir(vtk_dir) if f.endswith('.vtk')])

# Setup VTK renderer and window
renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
render_window.SetSize(640, 480)

# Setup window to image filter
window_to_image_filter = vtk.vtkWindowToImageFilter()
window_to_image_filter.SetInput(render_window)
window_to_image_filter.ReadFrontBufferOff()

# Setup video writer (OpenCV)
video_writer = None

for vtk_file in vtk_files:
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(os.path.join(vtk_dir, vtk_file))
    reader.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(reader.GetOutput())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    renderer.RemoveAllViewProps()
    renderer.AddActor(actor)
    renderer.ResetCamera()

    render_window.Render()
    window_to_image_filter.Modified()
    window_to_image_filter.Update()

    vtk_image = window_to_image_filter.GetOutput()
    width, height, _ = vtk_image.GetDimensions()

    vtk_array = vtk_image.GetPointData().GetScalars()
    components = vtk_array.GetNumberOfComponents()
    arr = numpy_support.vtk_to_numpy(vtk_array).reshape(height, width, components)
    frame = cv2.cvtColor(arr, cv2.COLOR_RGB2BGR)

    if video_writer is None:
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video_writer = cv2.VideoWriter(output_video, fourcc, frame_rate, (width, height))

    video_writer.write(frame)

if video_writer:
    video_writer.release()

print(f"Video saved as {output_video}")

import xarray as xr
import numpy as np
from vtk import vtkPoints, vtkPolyData, vtkCellArray, vtkPolyLine, vtkXMLPolyDataWriter, vtkIntArray

xsize = 3200   
ysize = 1600   
max_dxy = 1000
trail_length = 75

ds = xr.open_dataset('particle_dump.0000000.h5')
x = ds['x'].values
y = ds['y'].values
z = ds['z'].values
n_time, n_particles = x.shape

dx = np.diff(x, axis=0)
dy = np.diff(y, axis=0)

for t_end in range(1, n_time + 1):
    print(f'Parsing {t_end}/{n_time}')

    t_start = max(0, t_end-trail_length)
    n_steps = t_end-t_start

    points = vtkPoints()
    lines = vtkCellArray()
    pid_array = vtkIntArray()
    pid_array.SetName('particle_id')

    # Exclude cyclic BC crossings.
    for p in range(n_particles):
        if t_end > 1:
            dxi = dx[t_start:t_end-1, p]
            dyi = dy[t_start:t_end-1, p]
            if np.any(np.abs(dxi) > max_dxy) or np.any(np.abs(dyi) > max_dxy):
                continue

        polyline = vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(n_steps)
        for i, t in enumerate(range(t_start, t_end)):
            pid = points.InsertNextPoint(float(x[t, p]), float(y[t, p]), float(z[t, p]))
            polyline.GetPointIds().SetId(i, pid)
            pid_array.InsertNextValue(int(p))
        lines.InsertNextCell(polyline)

    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetPointData().AddArray(pid_array)
    polydata.GetPointData().SetActiveScalars('particle_id')

    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(f'particle_traj_t{t_end:04d}.vtp')
    writer.SetInputData(polydata)
    writer.Write()

# Generate x = 100 by y = 1 mesh for laglos shocktube experiments
import numpy as np

def main():
   nx_gridpoints = 15
   ny_gridpoints = 7
   xL = 0.
   xR = 7.
   yL = 0.
   yR = 3.
   x_arr = np.linspace(xL,xR,nx_gridpoints)
   y_arr = np.linspace(yL,yR,ny_gridpoints)
   home_dir = "/Users/madisonsheridan/Workspace/Laglos/"
   # filename = home_dir + "data/vortex-square-131044.mesh"
   filename = home_dir + "data/triple-point.mesh"
   f = open(filename, "w")

   # Prelimary information to write to mesh file
   f.write("MFEM mesh v1.0\n")
   f.write("\n")
   f.write("#\n")
   f.write("# MFEM Geometry Types (see mesh/geom.hpp):\n")
   f.write("#\n")
   f.write("# POINT       = 0\n")
   f.write("# SEGMENT     = 1\n")
   f.write("# TRIANGLE    = 2\n")
   f.write("# SQUARE      = 3\n")
   f.write("# TETRAHEDRON = 4\n")
   f.write("# CUBE        = 5\n")
   f.write("# PRISM       = 6\n")
   f.write("#\n")
   f.write("\n")
   f.write("dimension\n")
   f.write("2\n")
   f.write("\n")

   # ELEMENTS
   f.write("elements\n")
   # Num elements goes here
   f.write(str((nx_gridpoints-1)*(ny_gridpoints-1)) + "\n")
   # List of vertices for each element
   # el = 1
   for i in range(0, nx_gridpoints - 1):
      el_attr = 0
      for j in range(0, ny_gridpoints - 1):
         if (i < (nx_gridpoints - 1) / 7.):
            el_attr = 1
         elif (j%ny_gridpoints >= (ny_gridpoints - 1) / 2.):
            el_attr = 1
         else:
            el_attr = 2
         left = j + i*ny_gridpoints
         f.write("%d 3 %d %d %d %d\n" % (el_attr, left, left+ny_gridpoints, left+ny_gridpoints+1, left+1))
         # el += 1
   f.write("\n")

   # BOUNDARY
   f.write("boundary\n")

   bdry_bottom = 1
   bdry_right = 2
   bdry_top = 3
   bdry_left = 4
   
   # Num boundary faces
   f.write(str(2*(nx_gridpoints - 1) + 2 * (ny_gridpoints - 1)) + "\n")
   # TODO: Vertices connecting boundary edges

   #BOTTOM
   for i in range(0, (nx_gridpoints - 1) * ny_gridpoints, ny_gridpoints):
      f.write("%d %d %d %d\n" % (bdry_bottom, 1, i, i+ny_gridpoints))

   #RIGHT
   for i in range((nx_gridpoints-1) * ny_gridpoints, nx_gridpoints * ny_gridpoints - 1):
      f.write("%d 1 %d %d\n" % (bdry_right, i, i+1))

   #TOP
   for i in range(nx_gridpoints * ny_gridpoints - 1, ny_gridpoints - 1, -ny_gridpoints):
      f.write("%d %d %d %d\n" % (bdry_top, 1, i, i-ny_gridpoints))
   
   #LEFT
   for i in range(ny_gridpoints - 1, 0, -1):
      f.write("%d 1 %d %d\n" % (bdry_left, i, i-1))

   f.write("\n")

   # VERTICES
   f.write("vertices\n")
   f.write(str(nx_gridpoints * ny_gridpoints) + "\n")
   f.write("2\n")
   for i in range(0, nx_gridpoints):
      x = x_arr[i]
      for j in range(0, ny_gridpoints):
         y = y_arr[j]
         f.write("%.2f %.2f\n" % (x, y))

main()
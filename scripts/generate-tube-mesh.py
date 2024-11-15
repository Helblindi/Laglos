# Generate x = 100 by y = 1 mesh for laglos shocktube experiments
import numpy as np

def main():
   nx_gridpoints = 28
   ny_gridpoints = 5

   x_arr = np.linspace(-1.7,1.,nx_gridpoints)
   y_arr = np.linspace(0.,.4,ny_gridpoints)
   home_dir = "/Users/madisonsheridan/Workspace/Laglos/"
   filename = home_dir + "data/shocktube-test.mesh"
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
   # TODO: Num elements goes here
   f.write(str((nx_gridpoints-1) * (ny_gridpoints-1)) + "\n")
   # TODO: List of vertices for each element
   el = 1
   for i in range(0, (nx_gridpoints - 1) * ny_gridpoints, ny_gridpoints):
      for j in range(i, i + ny_gridpoints - 1):
         print("i: ", i)
         f.write("%d 3 %d %d %d %d\n" % (el, j, j+ny_gridpoints, j+ny_gridpoints+1, j+1))
         el += 1
   f.write("\n")

   # BOUNDARY
   f.write("boundary\n")
   bdry_bottom = 2
   bdry_right = 1
   bdry_top = 2
   bdry_left = 1
   
   # Num boundary faces
   f.write(str(2*(nx_gridpoints-1) + 2*(ny_gridpoints-1)) + "\n")
   # TODO: Vertices connecting boundary edges

   #BOTTOM
   for i in range(0, (nx_gridpoints - 1) * ny_gridpoints, ny_gridpoints):
      f.write("%d %d %d %d\n" % (bdry_bottom, 1, i, i+ny_gridpoints))

   #RIGHT
   for i in range((nx_gridpoints - 1) * ny_gridpoints, nx_gridpoints * ny_gridpoints - 1):
      f.write("%d %d %d %d\n" % (bdry_right, 1, i, i+1))

   #TOP
   for i in range(nx_gridpoints * ny_gridpoints - 1, ny_gridpoints, -ny_gridpoints):
      f.write("%d %d %d %d\n" % (bdry_top, 1, i, i-ny_gridpoints))

   #LEFT
   for i in range(ny_gridpoints - 1, 0, -1):
      f.write("%d %d %d %d\n" % (bdry_left, 1, i, i-1))

   f.write("\n")

   # VERTICES
   f.write("vertices\n")
   f.write(str(nx_gridpoints * ny_gridpoints) + "\n")
   f.write("2\n")
   for i in range(0, nx_gridpoints):
      for j in range(0, ny_gridpoints):
         # print("i: %.2f, j: %.2f" % (x_arr[i], y_arr[j]))
         f.write("%.6f %.6f\n" % (x_arr[i], y_arr[j]))

main()
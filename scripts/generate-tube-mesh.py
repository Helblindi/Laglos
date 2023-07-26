# Generate x = 100 by y = 1 mesh for laglos shocktube experiments
import numpy as np

def main():
   nx_gridpoints = 101
   ny_gridpoints = 2

   x_arr = np.linspace(0,1,nx_gridpoints)
   y_arr = np.linspace(0,1,ny_gridpoints)
   home_dir = "/Users/madisonsheridan/Workspace/Laglos/"
   filename = home_dir + "data/shocktube.mesh"
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
   f.write(str(int(nx_gridpoints * ny_gridpoints / 2) - 1) + "\n")
   # TODO: List of vertices for each element
   el = 1
   for i in range(0, (nx_gridpoints - 1) * ny_gridpoints, 2):
      print("i: ", i)
      f.write("%d 3 %d %d %d %d\n" % (el, i, i+2, i+3, i+1))
      el += 1
   f.write("\n")

   # BOUNDARY
   f.write("boundary\n")
   bdry_left = 0
   bdry_bottom = 1
   bdry_right = 2
   bdry_top = 3
   
   # Num boundary faces
   f.write(str(nx_gridpoints * ny_gridpoints) + "\n")
   # TODO: Vertices connecting boundary edges
   
   #LEFT
   bdr_el = 1
   f.write("%d 1 0 1\n" % bdr_el)
   bdr_el += 1

   #BOTTOM
   for i in range(0, (nx_gridpoints - 1) * ny_gridpoints, 2):
      f.write("%d %d %d %d\n" % (bdr_el, 1, i, i+2))
      bdr_el += 1

   #RIGHT
   f.write("%d 1 %d %d\n" % (bdr_el, (nx_gridpoints - 1) * ny_gridpoints, nx_gridpoints * ny_gridpoints - 1))
   bdr_el += 1

   #TOP
   for i in range(nx_gridpoints * ny_gridpoints - 1, 1, -2):
      f.write("%d %d %d %d\n" % (bdr_el, 1, i, i-2))
      bdr_el += 1

   f.write("\n")

   # VERTICES
   f.write("vertices\n")
   f.write(str(nx_gridpoints * ny_gridpoints) + "\n")
   f.write("2\n")
   for i in range(0, nx_gridpoints):
      for j in range(0, ny_gridpoints):
         # print("i: %.2f, j: %.2f" % (x_arr[i], y_arr[j]))
         f.write("%.2f %.2f\n" % (x_arr[i], y_arr[j]))

main()
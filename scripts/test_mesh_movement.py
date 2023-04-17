import numpy as np


def vel_exact(coords):
   # Parameters for velocity field
   a = .1
   b = 1.
   c = 1.
   d = 1.
   e = -0.1
   f = 1.
   # compute velocity components
   vx = a*coords[0] + b*coords[1] + c
   vy = d*coords[0] + e*coords[1] + f
   return np.array([vx, vy])

   
def identify_node_type(coords, index):
   # boundary nodes return 0
   if (coords[0] == 0 or coords[1] == 0 or coords[0] == 1 or coords[1] == 1):
      return 0
   # interior corner nodes return 1
   elif (coords[0] % 0.5 == 0 and coords[1] % 0.5 == 0):
      return 0
   # interior face nodes return 2
   elif ((coords[0] % 0.5 != 0 and coords[1] % 0.5 == 0) or
        (coords[0] % 0.5 == 0 and coords[1] % 0.5 != 0)):
      return index
   # cell center return 3
   else:
      return 0


def orthogonal(coords):
   x = -1 * coords[1]
   y = coords[0]
   return np.array([x,y])


def perp(coords):
   x = coords[1]
   y = -1 * coords[0]
   return np.array([x,y])


def compute_D(dt, V1, V2, n, A1, A2):
   t1 = np.subtract(V1, V2)
   nR = orthogonal(n)
   t2 = np.subtract(A2, A1)
   t2 = orthogonal(t2)
   D = dt * np.dot(t1, nR) + 2 * np.dot(n, t2)
   return D


def compute_bmn(Vf, nf, factor):
   val = np.dot(Vf, nf)
   val *= factor
   return val

def compute_c0(D, bmn, V1, V2, A1, A2, a3):
   A1R = orthogonal(A1)
   A2R = orthogonal(A2)
   a3R = orthogonal(a3)
   t1 = np.subtract(V2, V1)
   const1 = np.dot(V1, A1R) - np.dot(V2, A2R)
   const2 = np.dot(V1, A2R) - np.dot(V2, A1R)
   const3 = np.dot(t1, a3R)
   val = bmn + 0.5 * const1 + const2 / 6. + 2. * const3 / 3.
   val *= 3.
   val /= D
   return val


def compute_c1(dt, D, V1, V2, n, A1, A2):
   t1 = np.subtract(V2, V1)
   t2 = np.subtract(A2, A1)
   t2 = orthogonal(t2)
   n_perp = perp(n)
   val = dt * np.dot(t1, n) - 2 * np.dot(t2, n_perp)
   val /= D
   return val


def compute_V3perp(a3, a1new, a2new, a12new, c0, c1, dt, n):
   t1 = np.subtract(a2new, a1new)
   t2 = a3 - a12new + c0*dt*n
   t3 = c1 * n + perp(n)
   num = np.dot(t1, t2)
   denom = dt * np.dot(t1, t3)
   val = -num / denom
   return val


def main():
   # First create the mesh using a numpy array
   x = np.linspace(0,1,5)
   y = np.linspace(0,1,5)
   pos = np.zeros((np.size(x) * np.size(y), 5))
   faces = np.array([])
   face_normals = np.array([[0., 1.], [1.,0.], [1.,0.], [0.,1.]])
   face_adjacent_nodes = np.array([[2,12], [12, 10], [14, 12], [12, 22]])
   dt = 1.
   
   row = 0
   for i in range(0, np.size(x)):
      for j in range(0, np.size(y)): 
         print(row)
         pos[row,0] = x[i]
         pos[row,1] = y[j]
         # Get velocity
         pos[row,2:4] = vel_exact(pos[row,0:2])
         identifier = identify_node_type(pos[row, 0:2], row)
         pos[row, 4] = identifier
         if (identifier != 0):
            faces = np.append(faces, np.array([row]))
         
         # print(pos[index])
         row += 1
   
   for i in range(0, np.size(faces)):
      # compute face velocities as prescribed in the paper
      face_node = int(faces[i])
      print("-----------------------------------------")
      print("face: ", face_node)
      row = pos[face_node]
      face_x = row[0:2]
      face_v_exact = row[2:4]
      print("face_x: ", face_x)
      a1_index = face_adjacent_nodes[i,0]
      a2_index = face_adjacent_nodes[i,1]
      a1_x = np.array(pos[a1_index, 0:2])
      a1_v = np.array(pos[a1_index, 2:4])
      a2_x = np.array(pos[a2_index, 0:2])
      a2_v = np.array(pos[a2_index, 2:4])
      print("a1_x: ", a1_x)
      print("a2_x: ", a2_x)
      print("a1_v: ", a1_v)
      print("a2_v: ", a2_v)
      A1 = a1_x + 0.5*dt*a1_v
      A2 = a2_x + 0.5*dt*a2_v
      a1_full = a1_x + dt*a1_v
      a2_full = a2_x + dt*a2_v
      a12_full = 0.5 * (a1_full + a2_full)

      print("A1: ", A1)
      print("a1_full: ", a1_full)

      print("A2: ", A2)
      print("a2_full: ", a2_full)

      print("a12_full: ", a12_full)
      face_normal = face_normals[i]
      face_normal_R = orthogonal(face_normal)
      face_normal_perp = perp(face_normal)
      print("face_normal: ", face_normal)
      print("face_normal_R: ", face_normal_R)
      print("face_normal_perp: ", face_normal_perp)

      # Compute D
      D = compute_D(1., a1_v, a2_v, face_normal, A1, A2)
      print("D: ", D)

      # Compute half step normal and l2 norm of secant
      secant_half = np.subtract(A2, A1)
      factor = np.linalg.norm(secant_half)
      face_normal_half = orthogonal(secant_half)
      print("face normal half: ", face_normal_half)
      print("face norm factor: ", factor)

      # Compute bmn
      bmn = compute_bmn(face_v_exact, face_normal_half, 1.)
      print("bmn: ", bmn)

      # compute c0
      c0 = compute_c0(D, bmn, a1_v, a2_v, A1, A2, face_x)
      print("c0: ", c0)

      # compute c1
      c1 = compute_c1(dt, D, a1_v, a2_v, face_normal, A1, A2)
      print("c1: ", c1)

      # compute V3perp
      V3perp = compute_V3perp(face_x, a1_full, a2_full, a12_full, c0, c1, dt, face_normal)
      print("V3perp: ", V3perp)

      # compute V3n
      V3n = c1 * V3perp + c0
      print("V3n: ", V3n)

      # finally, compute the approximated velocity
      vel_approx = V3n * face_normal + V3perp * face_normal_perp
      print("vel_approx: ", vel_approx)

      print("face_v_exact: ", face_v_exact)


   # # display output of object
   # for row in pos:
   #    print(row)


main()
a = np.array([1,2])
print("norm of a: ", np.linalg.norm(a))
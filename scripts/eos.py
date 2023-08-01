
def main():
   a = 1.
   b = 1.
   rho = 4.9e-5
   gamma = 1.02
   pressure = 5.e-8
   sie = vdw_sie(gamma, rho, pressure, a, b)

   print("sie: ", sie)
   
    


def vdw_pressure(gamma, rho, sie, a, b):
   return (gamma - 1.) * (rho * sie + a * (rho**2)) / (1 - b * rho) - a * (rho**2)


def vdw_sie(gamma, rho, pressure, a, b):
   return ((1. / rho - b) * (pressure + a * (rho**2)) / (gamma - 1.)) - a * rho


main()
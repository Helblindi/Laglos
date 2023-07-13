import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
# ax = plt.axes(projection='3d')
# fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, projection='3d')
ax1 = fig.add_subplot(2,2,1,projection='3d')
ax2 = fig.add_subplot(2,2,2,projection='3d')
ax3 = fig.add_subplot(2,2,3,projection='3d')
ax4 = fig.add_subplot(2,2,4,projection='3d')


# Data for a three-dimensional line
# create x,y
nx = 100
ny = 100
x = np.linspace(0,1,nx)
y = np.linspace(0,1,ny)
xx, yy = np.meshgrid(x,y)

theta1 = 3./4. + (3./2.)*xx - (5./2.)*yy - (3./2.)*(xx**2 - yy**2)
theta2 = -1./4. - (1./2.)*xx + (3./2.)*yy + (3./2.)*(xx**2 - yy**2)
theta3 = -1./4. + (3./2.)*xx - (1./2.)*yy - (3./2.)*(xx**2 - yy**2)
theta4 = 3./4. - (5./2.)*xx + (3./2.)*yy + (3./2.)*(xx**2 - yy**2)

ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("$\Theta_1(\hat{x}, \hat{y})$")
ax1.plot_surface(xx, yy, theta1, alpha=0.5)

ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title("$\Theta_2(\hat{x}, \hat{y})$")
ax2.plot_surface(xx, yy, theta2, alpha=0.5)

ax3.set_xlabel("x")
ax3.set_ylabel("y")
ax3.set_title("$\Theta_3(\hat{x}, \hat{y})$")
ax3.plot_surface(xx, yy, theta3, alpha=0.5)

ax4.set_xlabel("x")
ax4.set_ylabel("y")
ax4.set_title("$\Theta_4(\hat{x}, \hat{y})$")
ax4.plot_surface(xx, yy, theta4, alpha=0.5)

plt.suptitle("Rotated Rannecher Turek Shape Functions")
plt.tight_layout()
plt.show()

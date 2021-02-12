import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

file =open("error_0.6_relglob.dat", "r")
lines = file.readlines()
x = np.zeros((1,1))
y = np.zeros((1,1))
z = np.zeros((1,1))

for line in lines:
    words= line.split(" ")

    if len(words) != 1:
        z = np.hstack([z, np.array([[float(words[2])]])])
        y = np.hstack([y, np.array([[float(words[1])]])])
        x = np.hstack([x, np.array([[float(words[0])]])])



# Make data.
X = np.arange(-1, 1, 0.5)
#print(X)
Y = np.arange(-1, 1, 0.5)
x, y = np.meshgrid(x, y)
X,Y = np.meshgrid(X,Y)
#z = np.meshgrid(z)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)
print(y)
print(x)
# Plot the surface.
#surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
                       #linewidth=0, antialiased=False)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
#ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

#plt.show()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(x, y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

#ax.scatter3D(x[1:, 0], y[1:, 0], z[1:, 0], cmap=cm.coolwarm, linewidth=0.05, antialiased=False)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
plt.title("Error")
plt.show()
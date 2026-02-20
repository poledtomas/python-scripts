import matplotlib.pyplot as plt
import numpy as np

dir = "/storage/brno12-cerit/home/poledto1/ising-hydro/smash-vhlle-hybrid/data/Hydro/"

col = "AuAu_7.7_2030"
data = np.loadtxt(dir + col + "/e_z0.dat")
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]
e = data[:, 3]


unique_x = np.unique(x)
unique_y = np.unique(y)
nx = len(unique_x)
ny = len(unique_y)

X = x.reshape((ny, nx))
Y = y.reshape((ny, nx))
E = e.reshape((ny, nx))

plt.figure(figsize=(8, 6))
plt.pcolormesh(X, Y, E, cmap="plasma")
plt.colorbar(label="e [GeV/fm^3]")
plt.xlabel("x[fm]")
plt.ylabel("y[fm]")
plt.title("Energy Density Distribution at z=0")
plt.savefig(f"{col}_heatmap.png", dpi=300)
plt.show()

plt.figure(figsize=(8, 6))
sc = plt.scatter(x, y, c=e, cmap="plasma", s=100)
plt.colorbar(sc, label="e [GeV/fm^3]")
plt.xlabel("x[fm]")
plt.ylabel("y[fm]")
plt.title("Energy Density Distribution at z=0")
plt.savefig(f"{col}_scatter.png", dpi=300)
plt.show()

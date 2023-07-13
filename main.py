import numpy as np
import matplotlib.pyplot as plt

def load_file(filename):
    with open(filename, "r") as txt_file:
        lines = txt_file.readlines()

    y = np.zeros(len(lines)-1)
    x = np.zeros(len(lines)-1)

    for i in range(0,len(lines)-1):
        x[i] = float(lines[i].split()[0])
        y[i] = float(lines[i].split()[1])

    return x, y

x,y = load_file("geo.txt")

xr,yr = load_file("out/r.txt")

x = x/1000
x = x - min(x) + 0.183
 
y = y/1000

print(np.pi*y[-1]**2)

plt.figure()
plt.axes().set_aspect('equal')
plt.plot(x,y)
plt.plot(xr,yr)
plt.show()


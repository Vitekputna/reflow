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

def load_comp(filename):
    with open(filename, "r") as txt_file:
        lines = txt_file.readlines()

    y1 = np.zeros(len(lines))
    y2 = np.zeros(len(lines))
    y3 = np.zeros(len(lines))

    x = np.zeros(len(lines))

    for i in range(0,len(lines)):
        x[i] = float(lines[i].split()[0])
        y1[i] = float(lines[i].split()[1])
        y2[i] = float(lines[i].split()[2])
        y3[i] = float(lines[i].split()[3])

    return x, y1,y2,y3

x,A = load_file("out/A.txt")
_,r = load_file("out/r.txt")
_,p = load_file("out/p.txt")
_,u = load_file("out/u.txt")
_,T = load_file("out/T.txt")
_,rho = load_file("out/W0.txt")
_,md = load_file("out/md.txt")
_,md_add = load_file("out/md_add.txt")
_,a = load_file("out/a.txt")
_,Y1,Y2,Y3 = load_comp("out/Y.txt")
_,H = load_file("out/H.txt")
_,H0 = load_file("out/H0.txt")

plt.figure(1)
plt.axes().set_aspect('equal')
plt.ylim([0,0.06])
plt.plot(x,r)
plt.title("Poloměr komory")
plt.xlabel("x[m]")
plt.ylabel("r[m]")
plt.grid()

plt.figure(2)
plt.plot(x,p/1e5)
plt.title("Tlak v komoře")
plt.xlabel("x[m]")
plt.ylabel("P[bar]")
plt.grid()

plt.figure(3)
plt.plot(x,u)
plt.title("Rychlost v komoře")
plt.xlabel("x[m]")
plt.ylabel(r"$u[ms^{-1}]$")
plt.grid()

plt.figure(4)
plt.plot(x,T)
plt.title("Teplota v komoře")
plt.xlabel("x[m]")
plt.ylabel("T[K]")
plt.grid()

plt.figure(5)
plt.plot(x,rho)
plt.title("Hustota v komoře")
plt.xlabel("x[m]")
plt.ylabel(r"$\rho[kg m^{-3}]$")
plt.grid()

plt.figure(6)
plt.plot(x,u/a)
plt.title("Machovo číslo v komoře")
plt.xlabel("x[m]")
plt.ylabel(r"$M[/]$")
plt.grid()

plt.figure(7)
plt.plot(x,Y1, label = "Produkty")
plt.plot(x,Y2, label = "Okysličovadlo")
plt.plot(x,Y3, label = "Palivo")
plt.title("Hmotnostní zlomky složek")
plt.xlabel("x[m]")
plt.ylabel(r"$Y[/]$")
plt.legend()
plt.grid()

plt.figure(8)
plt.plot(x,md)
plt.title("Hmotnostní tok")
plt.xlabel("x[m]")
plt.ylabel(r"$\dot{m}[kgs^{-1}]$")
plt.grid()

plt.figure(9)
plt.plot(x,H, label = "entalpie")
plt.plot(x,H0, label = "celková entalpie")
plt.title("Entalpie")
plt.xlabel("x[m]")
plt.ylabel(r"$H[Jm^{-3}]$")
plt.legend()
plt.grid()

plt.show()
import numpy as np
import matplotlib.pyplot as plt
import cmath as mth

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

def load(filename):
    with open(filename, "r") as txt_file:
        lines = txt_file.readlines()

    try:
        float(lines[0].split()[0])
    except:
        n_data = len(lines)-1
        start = 1
    else:
        n_data = len(lines)
        start = 0

    n_elements = len(lines[1].split())
    
    data = np.zeros((n_data,n_elements))

    for i in range(start,len(lines)):

        data[i-start,0] = float(lines[i].split()[0])
        
        for j in range(1,n_elements):
            data[i-start,j] = float(lines[i].split()[j])

    return data

def radius(X,N,rho):

    r = np.zeros((len(X),len(X[0])))

    N_sum = 0
    X_sum = 0
    try:
        for i in range(len(X)):

            for j in range(1,len(X[0])):

                N_sum += N[i,j]
                X_sum += abs(X[i,j])

                r[i,j] = (3*abs(X[i,j])/(4*mth.pi*rho*N[i,j]))**(1/3)
            
            
            r[i,0] = (3*X_sum/(4*mth.pi*rho*N_sum))**(1/3)
            N_sum = 0
            X_sum = 0
    except:
        return r
    else:
        return r

_,A = load_file("out/A.txt")
_,r = load_file("out/r.txt")
x,p = load_file("out/p.txt")
_,u = load_file("out/u.txt")
_,T = load_file("out/T.txt")
_,rho = load_file("out/W0.txt")
_,md = load_file("out/md.txt")
_,md_add = load_file("out/md_add.txt")
_,Y1,Y2,Y3 = load_comp("out/Y.txt")
_,H = load_file("out/H.txt")
_,H0 = load_file("out/H0.txt")
md_add = load("out/md_add.txt")
X = load("out/X.txt")
N = load("out/N.txt")
Ue = load("out/Ue.txt")
Te = load("out/Te.txt")
particles = load("out/particles.txt")

Ue = Ue[:,1:]
Te = Te[:,1:]

# plt.figure(1)
# plt.axes().set_aspect('equal')
# plt.ylim([0,0.06])
# plt.plot(x,r)
# plt.title("Poloměr komory")
# plt.xlabel("x[m]")
# plt.ylabel("r[m]")
# plt.grid()

plt.figure(2)
plt.plot(x,p/1e5)
plt.title("Tlak v komoře")
plt.xlabel("x[m]")
plt.ylabel("P[bar]")
plt.grid()

plt.figure(3)
plt.plot(x,u,'k',label = "Eulerian")
plt.plot(x,Ue,'b--',label = "Eulerian droplets")
plt.plot(particles[:,0],particles[:,2],'r-.',label = "Lagrangian droplets")
plt.title("Rychlost v komoře")
plt.xlabel("x[m]")
plt.ylabel(r"$u[ms^{-1}]$")
plt.legend()
plt.grid()

plt.figure(4)
plt.plot(x,T)
plt.plot(x,Te,'b--',label="Euler drop temperature")
plt.plot(particles[:,0],particles[:,3],'r',label = "Lagrangian")
plt.title("Teplota v komoře")
plt.xlabel("x[m]")
plt.ylabel("T[K]")
plt.grid()

# plt.figure(5)
# plt.plot(x,rho)
# plt.title("Hustota v komoře")
# plt.xlabel("x[m]")
# plt.ylabel(r"$\rho[kg m^{-3}]$")
# plt.grid()

# plt.figure(6)
# plt.plot(x,u/a)
# plt.title("Machovo číslo v komoře")
# plt.xlabel("x[m]")
# plt.ylabel(r"$M[/]$")
# plt.grid()

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

# plt.figure(9)
# plt.plot(x,H, label = "entalpie")
# plt.plot(x,H0, label = "celková entalpie")
# plt.title("Entalpie")
# plt.xlabel("x[m]")
# plt.ylabel(r"$H[Jm^{-3}]$")
# plt.legend()
# plt.grid()

plt.figure(10)
plt.plot(md_add[:,0],md_add[:,3])
plt.title("Hmotnostní tok z kapalné fáze")
plt.xlabel("x[m]")
plt.ylabel(r"$\dot{m}[kgs^{-1}]$")
plt.grid()

# plt.figure(11)
# plt.plot(X[:,0],X[:,1:])
# plt.title("Hmotnostní koncentrace kapalné fáze")
# plt.xlabel("x[m]")
# plt.ylabel(r"$X[kgm^{-3}]$")
# plt.grid()

# plt.figure(12)
# plt.plot(N[:,0],N[:,1:])
# plt.title("Početní koncentrace kapalné fáze")
# plt.xlabel("x[m]")
# plt.ylabel(r"$N[m^{-3}]$")
# plt.grid()

r = radius(X,N,700)
plt.figure(13)
plt.plot(x,r[:,1:])
plt.plot(x,r[:,0],'k--')
plt.plot(particles[:,0],particles[:,1],'r',label = "Lagrangian")
plt.title("Poloměr kapiček")
plt.xlabel("x[m]")
plt.ylabel(r"$r[m]$")
plt.grid()


plt.show()
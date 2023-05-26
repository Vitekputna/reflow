import numpy as np
import matplotlib.pyplot as plt
import cmath as mth
import sys

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
    
def plot_euler(path):
    _,r = load_file(path+"r.txt")
    x,p = load_file(path+"p.txt")
    _,u = load_file(path+"u.txt")
    _,T = load_file(path+"T.txt")
    _,rho = load_file(path+"W0.txt")
    _,md = load_file(path+"md.txt")
    _,md_add = load_file(path+"md_add.txt")
    _,Y1,Y2,Y3 = load_comp(path+"Y.txt")
    md_add = load(path+"md_add.txt")
    X = load(path+"X.txt")
    N = load(path+"N.txt")
    Ue = load(path+"Ue.txt")
    Te = load(path+"Te.txt")
    particles = load(path+"particles.txt")

    Ue = Ue[:,1:]
    Te = Te[:,1:]

    plt.figure(1)
    plt.plot(x,u,'k')
    plt.plot(x,Ue,'b--',label = "Eulerian droplets")
    plt.title("Rychlost v komoře")
    plt.xlabel("x[m]")
    plt.ylabel(r"$u[ms^{-1}]$")

    plt.figure(2)
    plt.plot(x,T,'k')
    plt.plot(x,Te,'b--',label="Eulerian droplets")
    plt.title("Teplota v komoře")
    plt.xlabel("x[m]")
    plt.ylabel("T[K]")
    plt.grid()

    r = radius(X,N,700)
    plt.figure(3)
    # plt.plot(x,r[:,1:],'')
    plt.plot(x,r[:,0],'b--',label = "Eulerian")
    plt.title("Poloměr kapiček")
    plt.xlabel("x[m]")
    plt.ylabel(r"$r[m]$")
    plt.grid()

    plt.figure(4)
    plt.plot(md_add[:,0],md_add[:,3],'b',label="Eulerian")
    plt.title("Hmotnostní tok z kapalné fáze")
    plt.xlabel("x[m]")
    plt.ylabel(r"$\dot{m}[kgs^{-1}]$")
    plt.grid()



def plot_lagrange(path):
    _,r = load_file(path+"r.txt")
    x,p = load_file(path+"p.txt")
    _,u = load_file(path+"u.txt")
    _,T = load_file(path+"T.txt")
    _,rho = load_file(path+"W0.txt")
    _,md = load_file(path+"md.txt")
    _,md_add = load_file(path+"md_add.txt")
    _,Y1,Y2,Y3 = load_comp(path+"Y.txt")
    md_add = load(path+"md_add.txt")
    X = load(path+"X.txt")
    N = load(path+"N.txt")
    Ue = load(path+"Ue.txt")
    Te = load(path+"Te.txt")
    particles = load(path+"particles.txt")

    Ue = Ue[:,1:]
    Te = Te[:,1:]

    plt.figure(1)
    plt.plot(x,u,'k')
    plt.plot(particles[:,0],particles[:,2],'r--',label = "Lagrangian droplets")
    plt.title("Rychlost v komoře")
    plt.xlabel("x[m]")
    plt.ylabel(r"$u[ms^{-1}]$")

    plt.figure(2)
    plt.plot(x,T,'k')
    plt.plot(particles[:,0],particles[:,3],'r--',label = "Lagrangian droplets")
    plt.title("Teplota v komoře")
    plt.xlabel("x[m]")
    plt.ylabel("T[K]")
    plt.grid()

    plt.figure(3)
    plt.plot(particles[:,0],particles[:,1],'r--',label = "Lagrangian")
    plt.title("Poloměr kapiček")
    plt.xlabel("x[m]")
    plt.ylabel(r"$r[m]$")
    plt.grid()

    plt.figure(4)
    plt.plot(md_add[:,0],md_add[:,3],'r',label="Lagrangian")
    plt.title("Hmotnostní tok z kapalné fáze")
    plt.xlabel("x[m]")
    plt.ylabel(r"$\dot{m}[kgs^{-1}]$")
    plt.grid()

plot_euler("tests/euler_all/")
plot_lagrange("tests/lagrange_all/")

plt.figure(1)
plt.grid()
plt.legend()

plt.figure(2)
plt.grid()
plt.legend()

plt.figure(3)
plt.grid()
plt.legend()

plt.show()
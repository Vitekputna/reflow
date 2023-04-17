import numpy as np
import matplotlib.pyplot as plt

def first_nonzero(data):
    
    for number in data:

        if number != 0:
            return number

def load_comp(filename):
    with open(filename, "r") as txt_file:
        lines = txt_file.readlines()

    n_elements = len(lines[1].split())
    n_data = len(lines)-1

    data = np.zeros((n_data,n_elements))

    for i in range(1,len(lines)):

        data[i-1,0] = float(lines[i].split()[0])
        
        for j in range(1,n_elements):
            data[i-1,j] = float(lines[i].split()[j])

    return data

data = load_comp("out/res.txt")

plt.yscale("log")
plt.plot(data[:,0],data[:,1]/first_nonzero(data[:,1]),label = "Density residual")
plt.plot(data[:,0],data[:,2]/first_nonzero(data[:,2]),label = "Momentum residual")
plt.plot(data[:,0],data[:,3]/first_nonzero(data[:,3]),label = "Energy residual")
plt.grid()
plt.legend()
plt.show()
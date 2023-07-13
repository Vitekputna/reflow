import numpy as np
import matplotlib.pyplot as plt
import time
import sys

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

path = sys.argv[1]

data = load_comp(path + "res.txt")

x = data[:,0]
y1 = data[:,1]/first_nonzero(data[:,1])
y2 = data[:,2]/first_nonzero(data[:,2])
y3 = data[:,3]/first_nonzero(data[:,3])


plt.ion()
figure, ax = plt.subplots()
line1, = ax.plot(x, y1,label = "Density residual")
line2, = ax.plot(x, y2,label = "Momentum residual")
line3, = ax.plot(x, y3,label = "Energy residual")


plt.yscale("log")
plt.grid()
plt.legend()

_,x_max = ax.get_xlim()
y_min, y_max = ax.get_ylim()

while True:
    data = load_comp(path+"res.txt")

    x = data[:,0]
    y1 = data[:,1]/first_nonzero(data[:,1])
    y2 = data[:,2]/first_nonzero(data[:,2])
    y3 = data[:,3]/first_nonzero(data[:,3])

    line1.set_xdata(x)
    line1.set_ydata(y1)

    line2.set_xdata(x)
    line2.set_ydata(y2)

    line3.set_xdata(x)
    line3.set_ydata(y3)

    if max(x) > 0.9*x_max:
        ax.set_xlim([min(x),1.3*max(x)]) 
        _,x_max = ax.get_xlim()

    if min([min(y1),min(y2),min(y3)]) < 2*y_min:
        ax.set_ylim([0.8*y_min,y_max]) 
        y_min, y_max = ax.get_ylim()    


    figure.canvas.draw()
    figure.canvas.flush_events()

    time.sleep(0.1)
import matplotlib.pyplot as plt
import numpy as np
import os
#calculate analytical result
def analytical(x):
    u = 0.5*x*x - 0.5*x
    return u

x_analytical = np.linspace(0,1,50);
u_analytical = []
for x in x_analytical:
    u_analytical.append(analytical(x))

# import numerical result
filedir = "./output.txt"

x_numerical = []
y_numerical = []

with open(filedir,"r+") as outputIo:
    lines = outputIo.readlines()
    for line in lines:
        infoVec = line.split(" ")
        x_numerical.append(float(infoVec[1]))
        y_numerical.append(float(infoVec[2]))
    #outputIo.close()

#  plot the result
fig, ax = plt.subplots(figsize=(7, 5), tight_layout=True)
ax.plot(x_analytical,u_analytical,label="analytical solution")
ax.scatter(x_numerical,y_numerical,label="numerical solution",color='red')
#ax.plot(x_numerical,y_numerical,label="numerical solution")
ax.grid()
ax.set_xlim(0, 1)
ax.set_ylim(-0.15, 0)
ax.set_xlabel("x",fontsize=15)
ax.set_ylabel("u",fontsize=15)
plt.legend()
plt.show()


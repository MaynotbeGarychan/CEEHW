import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
#calculate analytical result
def funcSum(x,y):
    retVal = 0
    for n in range(1,1):
        retVal += (1/(2*n-1))*math.sin((2*n-1)*math.pi*x/4)*math.pow(math.e,-((2*n-1)*math.pi*y)/4)
    return retVal

def analytical(x,y):
    u = (4/math.pi)*funcSum(x,y)
    return u

def testz(x,y):
    return x + y

fig = plt.figure()
ax = fig.gca(projection='3d')
# analytical data
x_analytical_set = np.linspace(0,1,100)
y_analytical_set = np.linspace(0,2,100)
z_analytical_set = []
for x in x_analytical_set:
    for y in y_analytical_set:
        z_analytical_set.append(analytical(x,y))

fig = plt.figure()
ax = plt.axes(projection='3d')
X,Y= np.meshgrid(x_analytical_set,y_analytical_set)
Z = np.ones([100,100])
for i in range(100):
    for j in range(100):
        Z[i][j] = analytical(x_analytical_set[i],y_analytical_set[j])
#ax.plot_wireframe(X,Y,Z)

#numerical
filedir = "./output.txt"

x_numerical = []
y_numerical = []
z_numerical = []

with open(filedir,"r+") as outputIo:
    lines = outputIo.readlines()
    for line in lines:
        infoVec = line.split(" ")
        x_numerical.append(float(infoVec[1]))
        y_numerical.append(float(infoVec[2]))
        z_numerical.append(float(infoVec[3]))
        ax.scatter(float(infoVec[1]),float(infoVec[2]),float(infoVec[3]),color='r')
    #outputIo.close()
ax.plot_wireframe(X,Y,Z)
ax.set_title('analytical vs numerical')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u')



plt.show()

"""
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
"""
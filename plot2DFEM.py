import math
import numpy as np

# =====================================
#    analytical
# =====================================
#calculate analytical result
def funcSum(x,y):
    retVal = 0
    for n in range(1,100):
        retVal += (1/(2*n-1))*math.sin((2*n-1)*math.pi*x/4)*math.pow(math.e,-((2*n-1)*math.pi*y)/4)
    return retVal

def analytical(x,y):
    u = (4/math.pi)*funcSum(x,y)
    return u

# analytical data
x_analytical_set = np.linspace(0,1,100)
y_analytical_set = np.linspace(0,2,100)
z_analytical_set = []
for x in x_analytical_set:
    for y in y_analytical_set:
        z_analytical_set.append(analytical(x,y))

# =====================================
#    numerical
# =====================================
filedir = "./output.txt"

nodeid_numerical = []
x_numerical = []
y_numerical = []
z_numerical = []

with open(filedir,"r+") as outputIo:
    lines = outputIo.readlines()
    for line in lines:
        infoVec = line.split(" ")
        nodeid_numerical.append(int(infoVec[0][1:]))
        x_numerical.append(float(infoVec[1]))
        y_numerical.append(float(infoVec[2]))
        z_numerical.append(float(infoVec[3]))
    #outputIo.close()

# =====================================
#    plot the result
# =====================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def heatmap(x_list,y_list,z_list):
    data = pd.DataFrame(data={'x':x_list, 'y':y_list, 'z':z_list})
    data = data.pivot(index='x', columns='y', values='z')
    sns.heatmap(data)

def threeDscatter(ax,x_list,y_list,z_list):
    for i in range(len(x_list)):
        ax.scatter(x_list[i], y_list[i], z_list[i], color='r')

fig = plt.figure()
ax = fig.gca(projection='3d')
ax = plt.axes(projection='3d')
threeDscatter(ax,x_numerical,y_numerical,z_numerical)
ax.set_title('numerical')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u')
plt.show()

heatmap(x_numerical,y_numerical,z_numerical)

plt.show()

"""
    output to excel
"""
# collect the value

import xlwt
workbook = xlwt.Workbook(encoding = 'ascii')
worksheet = workbook.add_sheet('My Worksheet')
worksheet.write(0, 0, 'node index')
worksheet.write(0, 1, 'analytical solution')
worksheet.write(0, 2, 'numerical solution')
for i in range(len(nodeid_numerical)):
    worksheet.write(i+1, 0, f'x{nodeid_numerical[i]}')
    worksheet.write(i + 1, 1, f'{analytical(x_numerical[i],y_numerical[i])}')
    print(analytical(x_numerical[i],y_numerical[i]))
    worksheet.write(i + 1, 2, f'{z_numerical[i]}')
workbook.save('formatting.xls')

import numpy as np
import matplotlib.pyplot as plt
import random
import xlwt
import pandas as pd
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

def main_makeNode():
    nodeidList = []
    nodexList = []
    nodeyList = []
    xmin = 0
    xmax = 1.0
    ymin = 0
    ymax = 1.0
    xNum = 8
    yNum = 8

    makeUniformNode(xmin,xmax,xNum,ymin,ymax,yNum,nodeidList,nodexList,nodeyList)
    #makeBoudaryNode(xmin, xmax, dx, ymin, ymax, dy, nodeidList, nodexList, nodeyList)
    #makeRondomNodeList(xmin, xmax, ymin, ymax, 16, nodeidList, nodexList, nodeyList)

    fig,ax=plt.subplots()
    ax.set_title("Node")
    ax.set_aspect(1)
    ax.scatter(nodexList,nodeyList,color='red')
    plt.show()

    txtDir = r'C:\Users\Gary\Desktop\CEE2021\reportOneSrc\report3NodeSeeds2.txt'
    outputNodeSeeds(txtDir,nodeidList,nodexList,nodeyList)

def main_plotElem():
    file_dir = r'C:\Users\Gary\Desktop\CEE2021\reportOneSrc\report3Mesh2.txt.txt'
    nodelist, elemlist = readMesh(file_dir)
    fig, ax = plt.subplots()
    ax.set_title("Mesh")
    ax.set_aspect(1)
    for node in nodelist:
        ax.scatter(node[1],node[2],color='r',s=50)
    for elem in elemlist:
        nodeid1 = elem[1]
        nodeid2 = elem[2]
        nodeid3 = elem[3]
        plotOneElem(ax,nodelist[nodeid1-1][1],nodelist[nodeid1-1][2],
                   nodelist[nodeid2-1][1],nodelist[nodeid2-1][2],
                   nodelist[nodeid3-1][1],nodelist[nodeid3-1][2])
    plt.show()

def main_addBoundaryConditions():
    file_dir = r'report3Mesh2.txt'
    nodeList,elemList = readMesh(file_dir)
    boud_dir = r'report3Mesh2Boud.txt'
    boudNum = 0
    with open(boud_dir,'w+') as txt:
        for node in nodeList:
            nodeId = node[0]
            nodeX = node[1]
            nodeY = node[2]
            if nodeX == 0:
                txt.write(f'boud {nodeId} 0\n')
                boudNum = boudNum + 1
            elif nodeY == 0:
                txt.write(f'boud {nodeId} 0\n')
                boudNum = boudNum + 1
            elif nodeX == 1:
                txt.write(f'boud {nodeId} {nodeY}\n')
                boudNum = boudNum + 1
            elif nodeY == 1:
                txt.write(f'boud {nodeId} {nodeX}\n')
                boudNum = boudNum + 1
    print(boudNum)

def main_plotResults():
    CGnodeidList,CGxList,CGyList, CGvalList = readResultTxt('./output2CG.txt')
    GaussnodeidList, GaussxList, GaussyList, GaussvalList = readResultTxt('output2Gauss.txt')

    #outputToExcel(GaussnodeidList,GaussxList,GaussyList,GaussvalList,'GaussTable')
    outputToExcel(CGnodeidList,CGxList,CGyList,CGvalList,'CGTable')

    #threeDScatterPlot(CGxList,CGyList,CGvalList,'Conjugate Gradient Method')
    #threeDScatterPlot(GaussxList, GaussyList, GaussvalList, 'Gauss Pivot Elimination Method')
    #analyticalValList = []
    #for i in range(len(GaussnodeidList)):
    #    analyticalValList.append(formulaAnalytical(GaussxList[i],GaussyList[i]))
    #threeDScatterPlot(GaussxList,GaussyList,analyticalValList,'Analytical')

def formulaAnalytical(x,y):
    val = x * y
    return val

def outputToExcel(nodeidList,xList,yList,uList,nameStr):
    workbook = xlwt.Workbook(encoding='ascii')
    worksheet = workbook.add_sheet('My Worksheet')
    worksheet.write(0, 0, 'Node index')
    worksheet.write(0, 1, 'Analytical solution')
    worksheet.write(0, 2, 'Numerical solution')
    for i in range(len(nodeidList)):
        worksheet.write(i + 1, 0, f'x{nodeidList[i]}')
        worksheet.write(i + 1, 1, f'{formulaAnalytical(xList[i],yList[i])}')
        worksheet.write(i + 1, 2, f'{uList[i]}')
    workbook.save(nameStr+'.xls')

def threeDScatterPlot(x_list,y_list,z_list,nameStr):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for i in range(len(x_list)):
        ax.scatter(x_list[i], y_list[i], z_list[i], color='blue')
    ax.set_title(nameStr)
    plt.show()
    #plt.savefig(nameStr+'.jpg')

def readResultTxt(result_dir):
    nodeidList = []
    xList = []
    yList = []
    valList = []
    with open(result_dir,'r') as txt:
        lines = txt.readlines()
        for line in lines:
            infoVec = line.split(" ")
            nodeidList.append(int(infoVec[0][1:]))
            xList.append(float(infoVec[1]))
            yList.append(float(infoVec[2]))
            valList.append(float(infoVec[3]))
    return nodeidList,xList,yList,valList

def makeUniformNode(xmin,xmax,xNum,ymin,ymax,yNum,nodeidList,nodexList,nodeyList):
    dx = (xmax-xmin)/xNum
    dy = (ymax-ymin)/yNum
    xlist = []
    ylist = []
    x = 0
    y = 0
    for i in range(xNum+1):
        xlist.append(x)
        x = x + dx
    for i in range(yNum+1):
        ylist.append(y)
        y = y + dy
    nodeid = 1
    for i in range(xNum+1):
        for j in range(yNum+1):
            nodeidList.append(nodeid)
            nodexList.append(xlist[i])
            nodeyList.append(ylist[j])
            nodeid = nodeid + 1

def plotOneElem(ax,x1,y1,x2,y2,x3,y3):
    ax.plot([x1, x2], [y1, y2], color='blue')
    ax.plot([x1, x3], [y1, y3], color='blue')
    ax.plot([x3, x2], [y3, y2], color='blue')

def readMesh(file_dir):
    nodeList = []
    elemList = []
    with open(file_dir,'r') as txt:
        lines = txt.readlines()
        for line in lines:
            infoVec = line.split()
            if infoVec[0] == 'mesh':
                nodeNum = int(infoVec[2])
                elemNum = int(infoVec[3])
            if infoVec[0] == 'node':
                nodeList.append([int(infoVec[1]),float(infoVec[2]),float(infoVec[3])])
            if infoVec[0] == 'elem':
                elemList.append([int(infoVec[1]),int(infoVec[2]),int(infoVec[3]),int(infoVec[4])])
    txt.close()
    return nodeList,elemList

def outputNodeSeeds(txtDir,nodeidList,nodexList,nodeyList):
    with open(txtDir,'w') as txt:
        for i in range(len(nodeidList)):
            txt.write(f"node {nodeidList[i]} {nodexList[i]} {nodeyList[i]}\n")
    txt.close()
    with open(txtDir, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(f'mesh 1 {len(nodeidList)} {0} {0}\n'+content)
    f.close()

def makeBoudaryNode(xmin,xmax,dx,ymin,ymax,dy,nodeidList,nodexList,nodeyList):
    xlist = []
    ylist = []
    x = xmin
    y = ymin
    xNum = int((xmax-xmin)/dx)+1
    yNum = int((ymax-ymin)/dy)+1
    for i in range(xNum):
        xlist.append(x)
        x = x + dx
    for i in range(yNum):
        ylist.append(y)
        y = y + dy
    totalNum = xNum *2 +(yNum-2)*2

    nodeid = 1
    for i in range(xNum):
        nodeidList.append(nodeid)
        nodexList.append(xlist[i])
        nodeyList.append(ymin)
        nodeid = nodeid +1
    for i in range(xNum):
        nodeidList.append(nodeid)
        nodexList.append(xlist[i])
        nodeyList.append(ymax)
        nodeid = nodeid +1
    for i in range(1,yNum-1):
        nodeidList.append(nodeid)
        nodexList.append(xmin)
        nodeyList.append(ylist[i])
        nodeid = nodeid + 1
    for i in range(1,yNum-1):
        nodeidList.append(nodeid)
        nodexList.append(xmax)
        nodeyList.append(ylist[i])
        nodeid = nodeid + 1

def makeRondomNodeList(xmin,xmax,ymin,ymax,num,nodeidList,nodexList,nodeyList):
    beginId = nodeidList[-1]+1
    for i in range(num):
        nodeidList.append(beginId)
        x,y = retRandomNodeCoor(xmin,xmax,ymin,ymax)
        nodexList.append(x)
        nodeyList.append(y)
        beginId = beginId + 1

def retRandomNodeCoor(xmin,xmax,ymin,ymax):
    #x = xmin + random.random() * (xmax - xmin)
    x = random.uniform(xmin,xmax)
    #y = ymin + random.random() * (ymax - ymin)
    y = random.uniform(ymin,ymax)
    return x,y


if __name__ == "__main__":
    #main_plotElem()
    #main_makeNode()
    #main_addBoundaryConditions()
    main_plotResults()
import math

mesh_dir = r'..\tmpInputOutput\report5Mesh.txt'
def main():
    nodeList, elemList = makeMesh(20, 1)
    numNode = len(nodeList)
    numElem = len(elemList)

    endTime = float(20)
    dt = float(0.4)
    dynaBoudList = makeDynamicBoudary(endTime, dt)
    statBoudList = makeStaticBoudary(nodeList)
    numDynaBoud = len(dynaBoudList)
    numStatBoud = len(statBoudList)
    numBoud = numDynaBoud + numStatBoud

    with open(mesh_dir,"w") as txt:
        txt.write(f'mesh 1 1 line {numNode} {numElem} {numDynaBoud}\n')
        for node in nodeList:
            txt.write(f'node {node[0]} {node[1]}\n')
        for elem in elemList:
            txt.write(f'elem {elem[0]} {elem[1]} {elem[2]}\n')
        txt.write(f'boudhead 1 {0} {numDynaBoud}\n')
        #for statBoud in statBoudList:
        #    txt.write(f'statboud {statBoud[0]} {statBoud[1]}\n')
        for dynaBoud in dynaBoudList:
            txt.write(f'dynaboud {dynaBoud[1]} 1\n')
            txt.write(f'dynaboud {dynaBoud[0]} {dynaBoud[2]}\n')
        txt.write(f'analysis 3\n')
        txt.write(f'analysis probName wave\n')
        txt.write(f'analysis solvName gpes\n')
        txt.write(f'analysis timeInte 1 0.0 {endTime} {dt} 0.25\n')
        txt.close()

def makeMesh(tolLen,ds):
    nodeList = []
    elemList= []

    numNode = int(tolLen/ds + 1)
    for i in range(numNode):
        nodeid = int(i+1)
        nodex = float(ds*i)
        nodeList.append([nodeid,nodex])

    numElem = int(tolLen/ds)
    for i in range(numElem):
        elemid = int(i + 1)
        leftNodeId = int(i + 1)
        rightNodeId = int(i + 2)
        elemList.append([elemid,leftNodeId,rightNodeId])

    return nodeList,elemList

def makeDynamicBoudary(timeLen,dt):
    nodeid = 1
    boudList=[]
    numTimeStep = int(timeLen/dt)

    for i in range(numTimeStep):
        currTime = i * dt
        boudList.append([nodeid,currTime,func(currTime)])

    return boudList

def func(t):
    alpha = 1/(2*math.pi)
    t0 = 10
    b = pow(alpha, 2) * pow(math.pi, 2) * pow(t - t0, 2)
    a = (1-2*b)*pow(math.e,-b)
    return a

def makeStaticBoudary(nodeList):
    boudList = []
    for i in nodeList:
        nodeId = int(i[0])
        val = float(0)
        boudList.append([nodeId,val])

    return boudList



if __name__ == "__main__":
    main()
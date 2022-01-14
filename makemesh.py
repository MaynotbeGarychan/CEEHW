import numpy as np

mesh_dir = r"report2Mesh3.txt"
"""
    specify the node
"""
xcoor_list=np.linspace(0,1,4)
ycoor_list=np.linspace(0,2,7)
#xcoor_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
#ycoor_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
xcoor_num = len(xcoor_list)
ycoor_num = len(ycoor_list)
"""
    make mesh
"""
with open(mesh_dir,"w+") as txt:
    # node
    node_num = len(xcoor_list)*len(ycoor_list)
    nodeid = 1
    for i in ycoor_list:
        for j in xcoor_list:
            txt.write(f"node {nodeid} {j} {i}\n")
            nodeid = nodeid + 1

    # elem
    elemid=1
    temp = 1
    nodeidMat = np.zeros((ycoor_num,xcoor_num))
    for i in range(ycoor_num):
        for j in range(xcoor_num):
            nodeidMat[i][j] = temp
            temp = temp +1
    for j in range(ycoor_num-1):
        for i in range(xcoor_num-1):
            txt.write(f"elem {elemid} {int(nodeidMat[j][i])} {int(nodeidMat[j][i+1])} {int(nodeidMat[j+1][i+1])}\n")
            elemid = elemid + 1
            txt.write(f"elem {elemid} {int(nodeidMat[j][i])} {int(nodeidMat[j+1][i + 1])} {int(nodeidMat[j+1][i])}\n")
            elemid = elemid + 1

    # boundary
    written_nodes_id = [] # specify a list that you would not duplicate boundary conditions
    boundaryid = 1
    for i in range(xcoor_num):
        txt.write(f"boun {int(nodeidMat[0][i])} 0\n")
        boundaryid = boundaryid+1
    for i in range(xcoor_num):
        txt.write(f"boun {int(nodeidMat[ycoor_num-1][i])} 1\n")
        boundaryid = boundaryid + 1

#    for i in range(ycoor_num):
#        if nodeidMat[i][0] not in written_nodes_id:
#            txt.write(f"boun {int(nodeidMat[i][0])} {pow(ycoor_list[i],2)}\n")
#            written_nodes_id.append(nodeidMat[i][0])
#            boundaryid = boundaryid + 1
#    for i in range(xcoor_num):
#        if nodeidMat[0][i] not in written_nodes_id:
#            txt.write(f"boun {int(nodeidMat[0][i])} {pow(xcoor_list[i],2)}\n")
#            written_nodes_id.append(nodeidMat[0][i])
#            boundaryid = boundaryid + 1
#    for i in range(ycoor_num):
#        if nodeidMat[i][1] not in written_nodes_id:
#            txt.write(f"boun {int(nodeidMat[i][1])} {1+pow(ycoor_list[i],2)}\n")
#            written_nodes_id.append(nodeidMat[i][1])
#            boundaryid = boundaryid + 1
#    for i in range(xcoor_num):
#        if nodeidMat[2][i] not in written_nodes_id:
#            txt.write(f"boun {int(nodeidMat[2][i])} {4+pow(xcoor_list[i],2)}\n")
#            written_nodes_id.append(nodeidMat[2][i])
#            boundaryid = boundaryid + 1
    txt.close()


    # write header info
    elem_num = elemid-1
    with open(mesh_dir, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(f'mesh 1 {nodeid-1} {elem_num} {boundaryid-1}\n' + content)
    f.close()
"""
import matplotlib.pyplot as plt
fig,ax=plt.subplots()
elem_node_list = np.zeros((elem_num,3))
tempid = 0
for j in range(ycoor_num - 1):
    for i in range(xcoor_num - 1):
        elem_node_list[tempid][0] = nodeidMat[j][i]
        elem_node_list[tempid][1] = nodeidMat[j][i+1]
        elem_node_list[tempid][2] = nodeidMat[j+1][i+1]
        tempid = tempid + 1
        elem_node_list[tempid][0] = nodeidMat[j][i]
        elem_node_list[tempid][1] = nodeidMat[j + 1][i + 1]
        elem_node_list[tempid][2] = nodeidMat[j + 1][i]
        tempid = tempid + 1


def searchXidYid(nodeid,mat,lenx,leny):
    for i in range(leny):
        for j in range (lenx):
            if nodeid == mat[i][j]:
                return i,j
    return 0,0

def pltOneElem(ax,x1,x2,x3,y1,y2,y3):
    ax.plot([x1,x2],[y1,y2],color='blue')
    ax.plot([x1, x3], [y1, y3],color='blue')
    ax.plot([x3, x2], [y3, y2],color='blue')

for i in range(elem_num):
    node1 = elem_node_list[i][0]
    node2 = elem_node_list[i][1]
    node3 = elem_node_list[i][2]
    x1 = xcoor_list[searchXidYid(node1,nodeidMat,xcoor_num,ycoor_num)[1]]
    x2 = xcoor_list[searchXidYid(node2, nodeidMat, xcoor_num, ycoor_num)[1]]
    x3 = xcoor_list[searchXidYid(node3, nodeidMat, xcoor_num, ycoor_num)[1]]
    y1 = ycoor_list[searchXidYid(node1,nodeidMat,xcoor_num,ycoor_num)[0]]
    y2 = ycoor_list[searchXidYid(node2, nodeidMat, xcoor_num, ycoor_num)[0]]
    y3 = ycoor_list[searchXidYid(node3, nodeidMat, xcoor_num, ycoor_num)[0]]
    pltOneElem(ax,x1, x2, x3, y1, y2, y3)

ax.set_title('mesh')
ax.set_aspect(1)
plt.show()
"""


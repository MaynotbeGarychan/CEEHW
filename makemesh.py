import numpy as np

mesh_dir = r"report2Mesh2.txt"

xcoor_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
ycoor_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]

#xcoor_list = [0.0,0.5,1.0]
#ycoor_list = [0.0,0.5,1.0,1.5,2.0]
xcoor_num = len(xcoor_list)
ycoor_num = len(ycoor_list)

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
    boundaryid = 1
    for i in range(xcoor_num):
        txt.write(f"boun {int(nodeidMat[0][i])} 0\n")
        boundaryid = boundaryid+1
    for i in range(xcoor_num):
        txt.write(f"boun {int(nodeidMat[ycoor_num-1][i])} 1\n")
        boundaryid = boundaryid + 1

    txt.close()

    # write header info
    with open(mesh_dir, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(f'mesh 1 {nodeid-1} {elemid-1} {boundaryid-1}\n' + content)
    f.close()
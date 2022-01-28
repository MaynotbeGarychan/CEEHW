# CEEHW
Author: CHEN Jiawei, 2nd year master student at the University of Tokyo

### Features
- PDE problem: 1D Wave, 2D Poisson, Newmark beta time integration scheme, Absorbing boundary condition(to be implemented).
- Matrix solver: Gauss pivot elimination algorithm, Conjugate gradient algorithm, matrix basic operation.
- Mesh generation: 2D Delaunay triangulation.
- Database: basic and simple data format, database, I/O for FEM.

### PDE problem:
- 2D Poisson:
![]()

- 1D Wave:
![]()



### I/O
Input format:

```buildoutcfg
mesh meshId nodeNum elemNum boudaryNum
node nodeId nodeXCoor nodeYCoor
elem elemId elemNode1Id elemNode2Id elemNode3Id
boud nodeId boundaryValue
```
Input example:
```buildoutcfg
mesh 1 6 6 4
node 1 0.000000 0.000000
node 2 1.000000 0.000000
node 3 1.000000 1.000000
elem 1 1 2 3
boud 1 0.5
```


### Matrix solver:



### 2D Delauney triangulation:



### Report One : 1DFEM with Gaussian Elimination

Here is the program structure:
![](https://github.com/MaynotbeGarychan/CEE_hw/blob/main/figures/structureOf1DFDM.JPG)

Problem:
![](https://github.com/MaynotbeGarychan/CEE_hw/blob/main/figures/Problem1.JPG)
![](https://github.com/MaynotbeGarychan/CEE_hw/blob/main/figures/mesh2OfProblem1.JPG)

Result:
![](https://github.com/MaynotbeGarychan/CEE_hw/blob/main/figures/FigureOf1DFDMProblem1.JPG)

### Report Two : 2DFEM with Gaussian Elimination
### Solving 2D Possion equation
![](https://github.com/MaynotbeGarychan/CEE_hw/blob/main/figures/Problem2JPG.JPG)

### Program structure
![](https://github.com/MaynotbeGarychan/CEE_hw/blob/main/figures/structureOf2DFEM.JPG)

### Performances of different meshes
![](https://github.com/MaynotbeGarychan/CEE_hw/blob/main/figures/Figure1Of2DFEMReport2.JPG)


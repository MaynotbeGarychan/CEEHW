# CEEHW
Author: CHEN Jiawei, 2nd year master student at the University of Tokyo

### Features
This is a repo to manage the HW of one of my attended courses.
Here are the features of this repo: 
- PDE problem: 1D Wave, 2D Poisson, Newmark beta time integration scheme, Absorbing boundary condition(to be implemented).
- Matrix solver: Gauss pivot elimination algorithm, Conjugate gradient algorithm, matrix basic operation.
- Mesh generation: 2D Delaunay triangulation.
- Database: basic and simple data format, database, I/O for FEM.

### PDE problem:
- The routine to solve the FEM problem:
![](.\figures\structureOf2DFEM.PNG)

- Case One - 2D Poisson:
![](.\figures\CaseOne2DPoisson.PNG)

- Case Two - 1D Wave:
![](.\figures\CaseTwo1DWave.PNG)

- Case Three - 1D Wave Dynamic, absorbing boundary condition with time integration scheme:
  (To be implemented)

- Input and Output format:
Input:
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
Output:
```buildoutcfg
nodeName nodeXCoor nodeYCoor nodeVal
```

Output example:
```buildoutcfg
x6 0.700000 0.200000 0.210770
```

### Matrix solver:
Gauss Pivot Elimination Solver and Conjugate Gradient Solver can be used to solve the FEM problem in this repo.
```buildoutcfg
int gaussianEliminationSolveMatrix(matrix* A, matrixInt* indexVec, matrix *result);
void conjugateSolveMatrix(const matrix systemMatrix, double tolerance, matrix* result);
```
Also, some basic matrix operations were built in this repo.
Here is the comparison of analytical results, numerical results from Gauss elimination solver and conjugate gradient 
solver for the 2D Poisson problem.
![](./figures/ComparisonAnalyticalNumerical.JPG)

### 2D Delaunay triangulation:
Delaunay triangulation is a common but popular mesh generation algorithm for triangle
element. 2D mesh generator for triangle element is supported in this repo now.
If you want to use it, please change the mode of the solution in the main.c by:
```buildoutcfg
int testMode = FEM_TEST_MODE;
```
Then it will go to the main function for mesh generator:
```buildoutcfg
int meshDelauneyTest(void)
```
- Input and Output:
To do the Delaunay triangulation, the node seeds have to be input, here is the format:
```buildoutcfg
mesh meshid NodeNum 0 0
node nodeid nodeXCoor nodeYCoor
```
```buildoutcfg
mesh 1 6 0 0
node 1 0 0
node 2 1.0 0
node 3 1.0 1.0
node 4 0 1.0
node 5 0.3 0.8
node 6 0.7 0.2
```
The output format of Delaunay triangulation is 
exactly, the mesh used for FEM simulation. Also, remember the node sequence
of each element has to be clockwise for FEM simulation, which has been done in 
this program:
```buildoutcfg
mesh meshId nodeNum elemNum 0
node nodeId nodeXCoor nodeYCoor
elem elemId node1Id node2Id node3Id
```
```buildoutcfg
mesh 1 6 6 4
node 1 0.000000 0.000000
node 2 1.000000 0.000000
node 3 1.000000 1.000000
node 4 0.000000 1.000000
node 5 0.300000 0.800000
node 6 0.700000 0.200000
elem 1 5 3 4
elem 2 5 4 1
elem 3 6 3 5
elem 4 6 5 1
elem 5 6 1 2
elem 6 6 2 3
```
- Case One - A very simple 2D mesh:
![](./figures/MeshCaseOne.JPG)

- Case Two - 2D mesh in 1x1 area:
![](./figures/MeshCaseTwo.JPG)


#include <stdio.h>
#include <malloc.h>
#include "1DFDM.h"
#include <string.h>

//#define ProblemOne
#define ProblemTwo

int main()
{
    //      reading mesh file     //
    // -------------------------  //
#ifdef ProblemOne
    struct meshInfo meshInfoDb;
    FILE *fileIo = fopen("problem1.txt","rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    char readType[4];

    fscanf(fileIo,"%s %d %d %d %d",readType,&(meshInfoDb.id),&(meshInfoDb.nodeNum),&(meshInfoDb.elementNum),&(meshInfoDb.boundaryNum));

    struct node nodeDb[meshInfoDb.nodeNum];
    struct element elementDb[meshInfoDb.elementNum];
    struct boundary boundaryDb[meshInfoDb.boundaryNum];
    
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[0].id),&(nodeDb[0].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[1].id),&(nodeDb[1].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[2].id),&(nodeDb[2].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[3].id),&(nodeDb[3].x));
    
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[0].id),&(elementDb[0].nodeId[0]),&(elementDb[0].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[1].id),&(elementDb[1].nodeId[0]),&(elementDb[1].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[2].id),&(elementDb[2].nodeId[0]),&(elementDb[2].nodeId[1]));

    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[0].nodeId),&(boundaryDb[0].value));
    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[1].nodeId),&(boundaryDb[1].value));
    fclose(fileIo);
#endif
#ifdef ProblemTwo
    struct meshInfo meshInfoDb;
    FILE *fileIo = fopen("problem2.txt","rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    char readType[4];

    fscanf(fileIo,"%s %d %d %d %d",readType,&(meshInfoDb.id),&(meshInfoDb.nodeNum),&(meshInfoDb.elementNum),&(meshInfoDb.boundaryNum));

    struct node nodeDb[meshInfoDb.nodeNum];
    struct element elementDb[meshInfoDb.elementNum];
    struct boundary boundaryDb[meshInfoDb.boundaryNum];
    
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[0].id),&(nodeDb[0].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[1].id),&(nodeDb[1].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[2].id),&(nodeDb[2].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[3].id),&(nodeDb[3].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[4].id),&(nodeDb[4].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[5].id),&(nodeDb[5].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[6].id),&(nodeDb[6].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[7].id),&(nodeDb[7].x));
    fscanf(fileIo,"%s %d %lf",readType,&(nodeDb[8].id),&(nodeDb[8].x));
    
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[0].id),&(elementDb[0].nodeId[0]),&(elementDb[0].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[1].id),&(elementDb[1].nodeId[0]),&(elementDb[1].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[2].id),&(elementDb[2].nodeId[0]),&(elementDb[2].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[3].id),&(elementDb[3].nodeId[0]),&(elementDb[3].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[4].id),&(elementDb[4].nodeId[0]),&(elementDb[4].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[5].id),&(elementDb[5].nodeId[0]),&(elementDb[5].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[6].id),&(elementDb[6].nodeId[0]),&(elementDb[6].nodeId[1]));
    fscanf(fileIo,"%s %d %d %d",readType,&(elementDb[7].id),&(elementDb[7].nodeId[0]),&(elementDb[7].nodeId[1]));

    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[0].nodeId),&(boundaryDb[0].value));
    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[1].nodeId),&(boundaryDb[1].value));
    fclose(fileIo);
#endif


    // assemble related matrix //
    // ----------------------- //
    // assemble loadvector from boundary conditions
    matrix loadVector;
    initilizeMatrix(&loadVector, meshInfoDb.nodeNum,1);
    assembleLoadVector(meshInfoDb,elementDb,nodeDb,&loadVector);

    // assemble element matrix (checked)
    matrix elemMatrix[meshInfoDb.elementNum];
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        initilizeMatrix(&(elemMatrix[i]),2,2);
        assembleElementStiffnessMatrix(elementDb[i],nodeDb,&(elemMatrix[i]));
        //printMatrix(&(elemMatrix[i]));
    }

    // assemble global matrix
    matrix globalMatrix;
    initilizeMatrix(&globalMatrix,meshInfoDb.nodeNum,meshInfoDb.nodeNum);
    assembleGlobalStiffnessMatrix(meshInfoDb,boundaryDb,elementDb,nodeDb,elemMatrix,&globalMatrix);

    // combine global matrix and loadvector to a linear system
    matrix linearSystem;
    newCombineMatrixCol(globalMatrix,loadVector,&linearSystem);


    // solve the linear system by Gauss Elimination //
    // -------------------------------------------- //
    // apply boundary condition to linear system
    applyBoundaryCondtion(&linearSystem,meshInfoDb,boundaryDb);
    
    // make a index vector to track the result index during swap
    matrixInt resultIndex;
    allocateMatrixInt(&resultIndex,linearSystem.numRow,1);
    for (int i = 0; i < resultIndex.numRow; i++)
    {
        resultIndex.mat[i][0] = nodeDb[i].id;
    }
    
    // solve the linear system with gauss elemination solver
    gaussianEliminationFDM(&linearSystem,&resultIndex);
    //printMatrix(&linearSystem);

    // print the result
    for (int i = 1; i < linearSystem.numRow-1; i++)
    {
        printf("x%d = %lf,\n",resultIndex.mat[i][0],linearSystem.mat[i][linearSystem.numCol-1]);
    }

    // Output the result //
    // ----------------- //
    FILE *OutputIo = fopen("output.txt","w");
    if (OutputIo == NULL)
    {
        return 0;
    }
    for (int i = 0; i < linearSystem.numRow; i++) // print the value of internal nodes
    {
        int nodeId = resultIndex.mat[i][0];
        if (linearSystem.mat[i][i] != 0)
        {
            fprintf(fileIo,"x%d %lf %lf\n",nodeId,nodeDb[nodeId-1].x,linearSystem.mat[i][linearSystem.numCol-1]);
        }
    }
    for (int i = 0; i < meshInfoDb.boundaryNum; i++) // print the value of boundary nodes
    {
        int nodeId =boundaryDb[i].nodeId;
        fprintf(fileIo,"x%d %lf %lf\n",nodeId,nodeDb[nodeId-1].x,boundaryDb[i].value);
    }
    
    fclose(OutputIo);

    return 1;
}

/**
 * @brief Matrix related functions
 * 
 */
// allocate the matrix
void allocateMatrix(matrix *T, int numRow, int numCol)
{
    T->mat = (double **)malloc(numRow*sizeof(double*));
    for (int i = 0; i < numRow; i++)
    {
        T->mat[i] = (double *)malloc(numCol*sizeof(double));
    }
    T->numCol = numCol;
    T->numRow = numRow;
}

void allocateMatrixInt(matrixInt *T, int numRow, int numCol)
{
    T->mat = (int **)malloc(numRow*sizeof(int*));
    for (int i = 0; i < numRow; i++)
    {
        T->mat[i] = (int *)malloc(numCol*sizeof(int));
    }
    T->numCol = numCol;
    T->numRow = numRow;
}

// initialize the value in matrix
void initilizeMatrix(matrix *T, int numRow, int numCol)
{
    allocateMatrix(T,numRow,numCol);
    for (int i = 0; i < numRow; i++)
    {
        for (int j = 0; j < numCol; j++)
        {
            T->mat[i][j] = 0.0;
        } 
    }
}

// free the matrix
void freeMatrix(matrix *T)
{
    for (int i = 0; i < T->numCol; i++)
    {
        free(T->mat[i]);
    }
    free(T->mat);
}

// print the matrix
void printMatrix(matrix *T)
{
    for (int i = 0; i < T->numRow; i++)
    {
        for (int j = 0; j < T->numCol; j++)
        {
            printf("%lf,",T->mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// fill the matrix for test
void fillMatrix33Test(matrix *T)
{
    T->mat[0][0] = 1;
    T->mat[0][1] = 1;
    T->mat[0][2] = 3;
    T->mat[1][0] = 6;
    T->mat[1][1] = 1;
    T->mat[1][2] = 2;
    T->mat[2][0] = 4;
    T->mat[2][1] = 0;
    T->mat[2][2] = 1;
}
void fillMatrix31Test(matrix *T)
{
    T->mat[0][0] = 2;
    T->mat[1][0] = 8;
    T->mat[2][0] = 3;
}

// Direct method
int gaussianEliminationFDM(matrix *A, matrixInt *indexVec)
/*begin
    Target: solve ax = b, A = { a | b }
    Step:
        1. Forward elimination
        2. Backward substitution
    Flag:   1: True  0:False
end*/
{
    // forward elimination
    if (!forwardElimintationPivot(A, indexVec))
    {
        return 0;
    }
    //printMatrix(A);
    // backward substitution
    if (!backwardSubtitution(A))
    {
        return 0;
    }
    //printMatrix(A);
    // rounde the diagonal component
    if (!roundDiagonalComponent(A))
    {
        return 0;
    }
    return 1;
}

// foward elimination with pivot
int forwardElimintationPivot(matrix *A, matrixInt *indexVec)
/*begin
*   Note: the all-zero row should be skipped in the future
end*/
{
    for (int i = 0; i < A->numRow; i++)
    {
        // if 0, need to be swap
        if (A->mat[i][i] == 0)
        {
            // find the maximum pivot
            int targetPos = i;
            for (int j = i+1; j < A->numRow-1; j++)
            {
                if (mathAbs(A->mat[j][i]) > mathAbs(A->mat[targetPos][i]))
                {
                    targetPos = j;
                }
            }
            swapRowMatrix(A,targetPos,i);
            swapRowMatrixInt(indexVec,targetPos,i);
        }
        // if non zero, we can eliminte it
        for (int j = i+1; j < A->numRow; j++)
        {
            if (A->mat[j][i] == 0)
            {
                continue;
            }
            double ratio = A->mat[j][i]/A->mat[i][i];
            A->mat[j][i] = 0;
            for (int k = i+1; k < A->numCol; k++)
            {
                A->mat[j][k] = A->mat[j][k] - ratio * A->mat[i][k];
            }
        }
    }
    return 1;
}



// transpose two line of matrix (by linePos)
int transposeMatrix(matrix *A, int rowPosOne, int rowPosTwo)
{
    if (rowPosOne >= A->numRow || rowPosTwo >= A->numRow )
    {
        return 0;
    }

    matrix tempRowOne;
    initilizeMatrix(&tempRowOne, 1, A->numCol);
    getRowOfMatrix(*A, rowPosOne, &tempRowOne);
    matrix tempRowTwo;
    initilizeMatrix(&tempRowTwo, 1, A->numCol);
    getRowOfMatrix(*A, rowPosTwo, &tempRowTwo);

    putRowOfMatrix(tempRowOne,rowPosTwo,A);
    putRowOfMatrix(tempRowTwo,rowPosOne,A);

    return 1;
}

// get any row of matrix
int getRowOfMatrix(matrix A, int rowPos, matrix *row)
/*
*   no need to init the row
*/
{
    if (rowPos >= A.numRow || row->numCol != A.numCol)
    {
        return 0;
    }
    allocateMatrix(row,1,A.numCol);
    for (int i = 0; i < A.numCol; i++)
    {
        row->mat[0][i] = A.mat[rowPos][i];
    }
    return 1;
}

int getColOfMatrix(matrix A, int colPos, matrix *col)
/*
*   please allocate, no need to init the col
*/
{
    if (colPos >= A.numCol || col->numRow != A.numRow)
    {
        return 0;
    }
    allocateMatrix(col,A.numRow,1);
    for (int i = 0; i < A.numCol; i++)
    {
        col->mat[i][0] = A.mat[colPos][i];
    }
    return 1;
}

// put any row of matrix
int putRowOfMatrix(matrix row, int rowPos, matrix *A)
{
    if (rowPos >= A->numRow || row.numCol != A->numCol)
    {
        return 0;
    }
    for (int i = 0; i < A->numCol; i++)
    {
        A->mat[rowPos][i] = row.mat[0][i];
    }
    return 1;
}

int inverseMatrix(matrix *A)
{
    matrix inverseA;
    initializeIdentityMatrix(&inverseA);
    matrix C;
    if (!newCombineMatrixCol(*A,inverseA,&C))
    {
        return 0;
    }
    forwardElimination(&C);
    backwardSubtitution(&C);
    roundDiagonalComponent(&C);
    getBlockOfMatrix(C,0,C.numRow-1,A->numCol,2*A->numCol-1,A);
    return 1;
}

void getBlockOfMatrix(matrix A, int beginRowPos, int endRowPos, int beginColPos ,int endColPos, matrix *block)
/*
*   Target: get a block from a matrix
*   Input: A, begin, end     Output: block (no need to init)
*/
{
    int blockNumRow = endRowPos-beginColPos+1;
    int blockNumCol = endColPos-beginColPos+1;
    allocateMatrix(block, blockNumRow,blockNumCol);
    for (int i = 0; i < blockNumRow; i++)
    {
        for (int j = 0; j < blockNumCol; j++)
        {
            block->mat[i][j] = A.mat[beginRowPos+i][beginColPos+j];
        }
    }
}

int roundDiagonalComponent(matrix *A)
{
    if (A->numRow>A->numCol)
    {
        return 0;
    }
    for (int i = 0; i < A->numRow; i++)
    {
        if (A->mat[i][i] == 1 || A->mat[i][i] == 0)
        {
            continue;
        }
        for (int j = A->numRow; j < A->numCol; j++)
        {
            A->mat[i][j] = A->mat[i][j]/A->mat[i][i];
        }
        A->mat[i][i] = 1.0;
    }
    return 1;
}

int newCombineMatrixCol(matrix A, matrix B, matrix *C)
/*
*   Target: C = { A | B }
*   Input: A, B     Output: C (no need to init)
*/
{
    if (A.numRow != B.numRow)
    {
        return 0;
    }
    // begin to combine
    allocateMatrix(C,A.numRow,A.numCol+B.numCol);
    for (int i = 0; i < A.numRow; i++)
    {
        for (int j = 0; j < A.numCol; j++)
        {
            C->mat[i][j] = A.mat[i][j];
        }
        for (int k = 0; k < B.numCol; k++)
        {
            C->mat[i][A.numCol+k] = B.mat[i][k];
        }
    }
    return 1;
}

int forwardElimination(matrix *A)
{
    if (A->numCol < A->numRow)
    {
        return 0;
    }
    for (int i = 0; i < A->numRow; i++)
    {
        for (int j = i+1; j < A->numRow; j++)
        {
            if (A->mat[j][i] == 0)
            {
                continue;
            }
            double ratio = A->mat[j][i]/A->mat[i][i];
            A->mat[j][i] = 0;
            for (int k = i+1; k < A->numCol; k++)
            {
                A->mat[j][k] = A->mat[j][k] - ratio * A->mat[i][k];
            }
        }
    }
    return 1;
}

int backwardSubtitution(matrix *A)
{
    if (A->numCol < A->numRow)
    {
        return 0;
    }
    for (int i = A->numRow - 1; i > 0; i--)
    {
        if (A->mat[i][i] == 0)
        {
            continue;
        }
        for (int j = i-1; j >= 0; j--)
        {
            if (A->mat[j][i] == 0)
            {
                continue;
            }
            double ratio = A->mat[j][i]/A->mat[i][i];
            A->mat[j][i] = 0.0;
            for (int k = A->numRow; k < A->numCol; k++)
            {
                A->mat[j][k] = A->mat[j][k] - ratio*A->mat[i][k];
            }
        }
    }
    return 1;
}

void swapRowMatrix(matrix *A,int rowOnePos,int rowTwoPos)
{
    double temp;
    for (int i = 0; i < A->numCol; i++)
    {
        temp = A->mat[rowOnePos][i];
        A->mat[rowOnePos][i] = A->mat[rowTwoPos][i];
        A->mat[rowTwoPos][i] = temp;
    }
}

void swapRowMatrixInt(matrixInt *A,int rowOnePos,int rowTwoPos)
{
    int temp;
    for (int i = 0; i < A->numCol; i++)
    {
        temp = A->mat[rowOnePos][i];
        A->mat[rowOnePos][i] = A->mat[rowTwoPos][i];
        A->mat[rowTwoPos][i] = temp;
    }
}

void initializeIdentityMatrix(matrix *A)
{
    for (int i = 0; i < A->numRow; i++)
    {
        for (int j = 0; j < A->numCol; j++)
        {
            if (i==j)
            {
                A->mat[i][j] = 1;
            }
            else
            {
                A->mat[i][j] = 0;
            }
        }
    }
}

int mutipleMatrix(matrix A, matrix B, matrix *C)
/*
*   Target: A * B = C(No need to init)
*/
{
    if (A.numCol != B.numRow)
    {
        return 0;     
    }
    allocateMatrix(C,A.numRow,B.numCol);
    for (int i = 0; i < C->numRow; i++)
    {
        matrix row;
        getRowOfMatrix(A,i,&row);
        for (int j = 0; j < C->numCol; j++)
        {
            matrix col;
            getColOfMatrix(B,j,&col);
            C->mat[i][j] = dotProduct(row,col);
        }
    }
    return 1;
}

double dotProduct(matrix a, matrix b)
{
    if (a.numCol != b.numRow)
    {
        return 0;
    }
    double val = 0;
    for (int i = 0; i < a.numCol; i++)
    {
        val += a.mat[0][i]*b.mat[i][0];
    }
    return val;
}

int subtractMatrix(matrix A, matrix B, matrix *C)
{
    if (A.numRow != B.numRow || A.numCol != B.numCol)
    {
        return 0;
    }
    allocateMatrix(C,A.numRow,A.numCol);
    for (int i = 0; i < C->numRow; i++)
    {
        for (int j = 0; j < C->numCol; j++)
        {
            C->mat[i][j] = A.mat[i][j] - B.mat[i][j];
        }
    }
    return 1;
}

/**
 * @brief math related functions
 * 
 */
double mathAbs(double a)
{
    if (a>=0)
    {
        return a;
    }
    else
    {
        return -a;
    }
}
/*
int allocateArrayInt(int *array, int length)
{
    array = (int **)malloc(length*sizeof(int*));
    for (int i = 0; i < length; i++)
    {
        array[i] = (int *)malloc(length*sizeof(int));
    }
}
*/

/**
 * @brief Pre post process related funtions
 * 
 */

// read header of txt
int readMeshHeader(FILE *fileIo, struct meshInfo *meshInfoDb)
/*begin
    Target: Read the basic information of meshes
    Flag:   1: True  0:False
end*/
{
    printf("Begin to read the input file\n");

    // meshinfo
    char readType[4];
    fscanf(fileIo,"%s %d %d %d",readType,&(meshInfoDb->nodeNum),&(meshInfoDb->elementNum),&(meshInfoDb->boundaryNum));
    
    return 1;
}


// read all the mesh information
int readMeshFile(const char* fileName, struct meshInfo *meshInfoDb, struct node *nodeDb)
{
    FILE *fileIo = fopen(fileName,"rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    printf("Begin to read the mesh file\n");

    // meshinfo
    char readType[4];
    fscanf(fileIo,"%s %d %d %d",readType,&(meshInfoDb->nodeNum),&(meshInfoDb->elementNum),&(meshInfoDb->boundaryNum));

    // nodeinfo
    for (int i = 0; i < meshInfoDb->nodeNum; i++,nodeDb++)
    {
        fscanf(fileIo,"%s %d %lf",readType,&(nodeDb->id),&(nodeDb->x));
    }

    return 1;
}

/**
 * @brief Solver related function
 * 
 */

// Assemble load vector (checked)
int assembleLoadVector(struct meshInfo meshInfoDb, struct element elementDb[],
         struct node nodeDb[], matrix *loadVector)
/*begin
*   Target: assmble loadvector
*/
{
    initilizeMatrix(loadVector, meshInfoDb.nodeNum, 1);
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        double temp = 0.5*(nodeDb[elementDb[i].nodeId[0]-1].x - nodeDb[elementDb[i].nodeId[1]-1].x);
        loadVector->mat[elementDb[i].nodeId[0]-1][0] += temp;
        loadVector->mat[elementDb[i].nodeId[1]-1][0] += temp;
    }
    //printMatrix(loadVector);
    return 1;
}


// Assemble element stiffness matrix
int assembleElementStiffnessMatrix(struct element elementDb,struct node nodeDb[],matrix *elemMatrix)
{   
    // calculate the basic component
    double xleft = nodeDb[elementDb.nodeId[0]-1].x;
    double xright = nodeDb[elementDb.nodeId[1]-1].x;
    double k = 1/(xright-xleft);

    // diagonal
    elemMatrix->mat[0][0] = k;
    elemMatrix->mat[1][1] = k;
    // non-diagonal
    elemMatrix->mat[0][1] = -k;
    elemMatrix->mat[1][0] = -k;

    return 1;
}

// Assemble global stiffness matrix (checked)
int assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb,struct boundary boundaryDb[],struct element elementDb[],
                struct node nodeDb[], matrix elemMatrix[],matrix *globalMatrix)
{
    // from element stiffness matrix
    for (int elemId = 0; elemId < meshInfoDb.elementNum; elemId++)
    {
        int leftNodePos = elementDb[elemId].nodeId[0]-1;
        int rightNodePos = elementDb[elemId].nodeId[1]-1;
        globalMatrix->mat[leftNodePos][leftNodePos] += elemMatrix[elemId].mat[0][0];
        globalMatrix->mat[rightNodePos][rightNodePos] += elemMatrix[elemId].mat[1][1];
        globalMatrix->mat[leftNodePos][rightNodePos] += elemMatrix[elemId].mat[0][1];
        globalMatrix->mat[rightNodePos][leftNodePos] += elemMatrix[elemId].mat[1][0];
    }
    //printMatrix(globalMatrix);
    // from boundary conditions
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        int boundaryNodeId = boundaryDb[i].nodeId;
        int linkNodeId = findSameElementNode(meshInfoDb,elementDb,boundaryNodeId);
        double temp = 1/mathAbs(nodeDb[boundaryNodeId-1].x-nodeDb[linkNodeId-1].x);
        globalMatrix->mat[boundaryNodeId-1][boundaryNodeId-1] += -temp;
        globalMatrix->mat[boundaryNodeId-1][linkNodeId-1] += temp;
    }
    //printMatrix(globalMatrix);
    return 1;
}

int findSameElementNode(struct meshInfo meshInfoDb, struct element elementDb[],int nodeId)
{
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        int nodeList[2];
        for (int j = 0; j < 2; j++)
        {
            nodeList[j] = elementDb[i].nodeId[j];
        }
        if (nodeId == nodeList[0])
        {
            return nodeList[1];
        }
        else if (nodeId == nodeList[1])
        {
            return nodeList[0];
        }
    }
    return 0;
}

int solveFEMSystem(matrix LHSMatrix, matrix RHSVector, struct meshInfo meshInfoDb, struct boundary boundaryDb[],
                struct node nodeDb[],matrix* result)
{
    // reorder the nodelist
    matrixInt reorderList;
    reorderNodeList(meshInfoDb,boundaryDb,nodeDb,&reorderList);

    // reorder the global matrix and load vector
    matrix reorderLHSMatrix;
    allocateMatrix(&reorderLHSMatrix,LHSMatrix.numRow,LHSMatrix.numCol);
    matrix reorderRHSVector;
    allocateMatrix(&reorderRHSVector,RHSVector.numRow,RHSVector.numCol);
    for (int i = 0; i < meshInfoDb.nodeNum; i++)
    {
        reorderRHSVector.mat[i][0] = RHSVector.mat[reorderList.mat[0][i]][0];
        for (int j = 0; j < LHSMatrix.numCol; j++)
        {
            reorderLHSMatrix.mat[i][j] = LHSMatrix.mat[reorderList.mat[0][i]][reorderList.mat[0][j]];
        }
    }

    // cut from matrix
    matrix Kaa;
    matrix Kab;
    matrix fa;
    matrix ub;
    // Kaa
    int numInternalNode = meshInfoDb.nodeNum - meshInfoDb.boundaryNum;
    getBlockOfMatrix(reorderLHSMatrix,0,numInternalNode-1,0,numInternalNode-1,&Kaa);
    // Kab
    getBlockOfMatrix(reorderLHSMatrix,0,numInternalNode,numInternalNode,meshInfoDb.nodeNum-1,&Kab);
    // fa
    getBlockOfMatrix(reorderRHSVector,0,numInternalNode-1,0,0,&fa);
    // ub
    allocateMatrix(&ub,meshInfoDb.boundaryNum,1);
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        ub.mat[i][0] = boundaryDb[i].value;
    }
    
    // calculation
    inverseMatrix(&Kaa);
    matrix Kabub;
    mutipleMatrix(Kab,ub,&Kabub);
    matrix faKabub;
    subtractMatrix(fa,Kabub,&faKabub);
    mutipleMatrix(Kaa,faKabub,result);

    return 1;
}

void reorderNodeList(struct meshInfo meshInfoDb, struct boundary boundaryDb[],
        struct node nodeDb[], matrixInt *reorderList)
{
    allocateMatrixInt(reorderList,1,meshInfoDb.nodeNum);

    int boundaryNodeList[meshInfoDb.boundaryNum];
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        boundaryNodeList[i] = boundaryDb[i].nodeId;
    }

    //int reorderIntNodeList[meshInfoDb.nodeNum-meshInfoDb.boundaryNum];
    int reorderBouNodeList[meshInfoDb.boundaryNum];

    int internalIndex = 0;
    int boundaryIndex = 0;
    for (int i = 0; i < meshInfoDb.nodeNum; i++)
    {
        int state = 0;
        for (int j = 0; j < meshInfoDb.boundaryNum; j++)
        {
            if (nodeDb[i].id == boundaryNodeList[j])
            {
                state = 1;
                break;
            }
            state = 0;
        }
        if (state == 1) // for boundary ub
        {
            reorderBouNodeList[boundaryIndex] = i;
            boundaryIndex++;
        }
        else // for internal Kaa, Kab, fa
        {
            reorderList->mat[0][i] = i;
            internalIndex++;
        }
    }

    for (int i = internalIndex; i < meshInfoDb.nodeNum; i++, reorderList++)
    {
        reorderList->mat[0][i] = reorderBouNodeList[i-internalIndex];
    }
}

void applyBoundaryCondtion(matrix *linearSystem, struct meshInfo meshInfoDb,struct boundary boundaryDb[])
{
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        double val = boundaryDb[i].value;
        int pos = boundaryDb[i].nodeId-1;
        int lvPos = linearSystem->numCol-1;
        for (int j = 0; j < linearSystem->numRow; j++)
        {
            linearSystem->mat[j][lvPos] -= linearSystem->mat[j][pos]*val;
            linearSystem->mat[j][pos] = 0.0;
        }
    }
}
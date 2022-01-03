#include <stdio.h>
#include <malloc.h>
#include "1DFDM.h"
#include <string.h>

int main()
{
    // reading mesh file //
    // ----------------  //
    struct meshInfo meshInfoDb;
    FILE *fileIo = fopen("problem1.txt","rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    char readType[4];
    fscanf(fileIo,"%s %d %d %d",readType,&(meshInfoDb.nodeNum),&(meshInfoDb.elementNum),&(meshInfoDb.boundaryNum));
    
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
    // end read mesh file

    // assemble related matrix //
    // ----------------------- //
    // assemble loadvector from boundary conditions
    meshInfoDb.nodeNum = 4;
    matrix loadVector;
    initilizeMatrix(&loadVector, meshInfoDb.nodeNum,1);
    assembleLoadVector(meshInfoDb,elementDb,nodeDb,&loadVector);

    // assemble element matrix
    matrix elemMatrix[meshInfoDb.elementNum];
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        initilizeMatrix(&(elemMatrix[i]),2,2);
        assembleElementStiffnessMatrix(elementDb[i],nodeDb,&(elemMatrix[i]));
    }

    // assemble global matrix
    matrix globalMatrix;
    initilizeMatrix(&globalMatrix,meshInfoDb.nodeNum,meshInfoDb.nodeNum);
    assembleGlobalStiffnessMatrix(meshInfoDb,boundaryDb,elementDb,nodeDb,elemMatrix,&globalMatrix);

    // combine loadvector and global stiffness matrix
    



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

// initialize the value in matrix
void initilizeMatrix(matrix *T, int numRow, int numCol)
{
    allocateMatrix(T,numRow,numCol);
    for (int i = 0; i < numRow; i++)
    {
        for (int j = 0; j < numCol; j++)
        {
            T->mat[i][j] = 1.0;
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
int gaussianElimination(matrix *A,matrix *b)
/*begin
    Target: solve Ax = b
    Step:
        1. Forward elimination
        2. Backward substitution
    Flag:   1: True  0:False
end*/
{
    if (A->numRow != b->numRow)
    {
        return 0;
    }
    // forward elimination
    if (!gaussForwardElimination(A,b))
    {
        return 0;
    }
    // backward substitution
    if (!gaussbackwardSubstitution(A,b))
    {
        return 0;
    }
    return 1;
}

// Forward elimination
int gaussForwardElimination(matrix *A, matrix *b)
{
    if(A->numRow != b->numRow)
    {
        return 0;
    }
    if (b->numCol != 1)
    {
        // current just support load vector
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
            b->mat[j][0] = b->mat[j][0] - ratio * b->mat[i][0];
        }
    }
    return 1;
}

// backward substitution
int gaussBackwardSubstitution(matrix *A, matrix* b)
{
    if(A->numRow != b->numRow)
    {
        return 0;
    }
    if (b->numCol != 1)
    {
        // current just support load vector
        return 0;
    }
    for (int i = A->numRow - 1; i > 0; i--)
    {
        for (int j = i-1; j >= 0; j--)
        {
            if (A->mat[j][i] == 0)
            {
                continue;
            }
            double ratio = A->mat[j][i]/A->mat[i][i];
            A->mat[j][i] = 0.0;
            b->mat[j][0] = b->mat[j][0]-ratio*b->mat[i][0];
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
{
    if (rowPos >= A.numRow || row->numCol != A.numCol)
    {
        return 0;
    }
    for (int i = 0; i < A.numCol; i++)
    {
        row->mat[0][i] = A.mat[rowPos][i];
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

void roundDiagonalComponent(matrix *A)
{
    for (int i = 0; i < A->numRow; i++)
    {
        for (int j = A->numRow; j < A->numCol; j++)
        {
            A->mat[i][j] = A->mat[i][j]/A->mat[i][j];
        }
    }
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

void forwardElimination(matrix *A)
{
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
}

void backwardSubtitution(matrix *A)
{
    for (int i = A->numRow - 1; i > 0; i--)
    {
        for (int j = i-1; j >= 0; j--)
        {
            if (A->mat[j][i] == 0)
            {
                continue;
            }
            A->mat[j][i] = 0.0;
            double ratio = A->mat[j][i]/A->mat[i][i];
            for (int k = A->numRow; k < A->numCol; k++)
            {
                A->mat[j][k] = A->mat[j][k] - ratio*A->mat[i][k];
            }
        }
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

/**
 * @brief Pre post process related funtions
 * 
 */

// read header of txt
int readMeshHeader(const char* fileName, struct meshInfo *meshInfoDb)
/*begin
    Target: Read the basic information of meshes
    Flag:   1: True  0:False
end*/
{
    FILE *fileIo = fopen(fileName,"rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    printf("Begin to read the input file\n");

    // meshinfo
    char readType[4];
    fscanf(fileIo,"%s %d %d %d",readType,&(meshInfoDb->nodeNum),&(meshInfoDb->elementNum),&(meshInfoDb->boundaryNum));
    
    fclose(fileIo);
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

// Assemble load vector
int assembleLoadVector(struct meshInfo meshInfoDb, struct element elementDb[],
         struct node nodeDb[], matrix *loadVector)
/*begin
*   Target: assmble loadvector without boundary conditions
*/
{
    initilizeMatrix(loadVector, meshInfoDb.nodeNum, 1);
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        matrix *elemLoadVector;
        initilizeMatrix(loadVector, 2, 1);
        double temp = 0.5*(nodeDb[elementDb[i].nodeId[0]-1].x - nodeDb[elementDb[i].nodeId[1]-1].x);
        loadVector->mat[elementDb[i].nodeId[0]-1][0] += temp;
        loadVector->mat[elementDb[i].nodeId[1]-1][0] += -temp;
    }
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
    elemMatrix->mat[1][1] = -k;

    return 1;
}

// Assemble global stiffness matrix
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
    
    // from boundary conditions
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        int boundaryNodeId = boundaryDb[i].nodeId;
        int linkNodeId = findSameElementNode(meshInfoDb,elementDb,boundaryNodeId);
        double temp = 1/abs(nodeDb[boundaryNodeId-1].x-nodeDb[linkNodeId-1].x);
        globalMatrix->mat[boundaryNodeId-1][boundaryNodeId-1] = temp;
        globalMatrix->mat[boundaryNodeId-1][linkNodeId-1] = -temp;
    }
    return 1;
}

int findSameElementNode(struct meshInfo meshInfoDb, struct element elementDb[],int nodeId)
{
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        int nodeList[2] = elementDb[i].nodeId;
        if (nodeId == nodeList[0])
        {
            return nodeList[1];
        }
        else if (nodeId == nodeList[1])
        {
            return nodeList[0];
        }
        else
        {
            return 0;
        }
    }
}

int solveSystem(matrix LHSMatrix, matrix RHSVector, struct meshInfo meshInfoDb, struct boundary boundaryDb[],
                struct node nodeDb[],matrix* result)
{
    // reorder the nodelist
    int reorderList[meshInfoDb.nodeNum];
    reorderNodeList(meshInfoDb,boundaryDb,nodeDb,&reorderList);

    // reorder the global matrix and load vector
    matrix reorderLHSMatrix;
    allocateMatrix(&reorderLHSMatrix,LHSMatrix.numRow,LHSMatrix.numCol);
    matrix reorderRHSVector;
    allocateMatrix(&reorderRHSVector,RHSVector.numRow,RHSVector.numCol);
    for (int i = 0; i < meshInfoDb.nodeNum; i++)
    {
        reorderRHSVector.mat[i][0] = RHSVector.mat[reorderList[i]][0];
        for (int j = 0; j < LHSMatrix.numCol; j++)
        {
            reorderLHSMatrix.mat[i][j] = LHSMatrix.mat[reorderList[i]][reorderList[j]];
        }
    }

    matrix Kaa;
    matrix Kab;
    matrix fa;
    matrix ub;


}

void reorderNodeList(struct meshInfo meshInfoDb, struct boundary boundaryDb[],
        struct node nodeDb[], int *reorderList[])
{
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
            reorderList[internalIndex] = i;
            internalIndex++;
        }
    }

    for (int i = internalIndex; i < meshInfoDb.nodeNum; i++)
    {
        reorderList[i] = reorderBouNodeList[i-internalIndex];
    }
}
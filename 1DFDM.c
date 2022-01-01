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
    //initilizeMatrix(&loadVector, meshInfoDb.nodeNum,1);
    initilizeMatrix(&loadVector, meshInfoDb.nodeNum,1);
    //assembleLoadVector(meshInfoDb,boundaryDb,&loadVector);(wrong)

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
    assembleGlobalStiffnessMatrix(meshInfoDb,elementDb,elemMatrix,&globalMatrix);
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
    if (!forwardElimination(A,b))
    {
        return 0;
    }
    // backward substitution
    if (!backwardSubstitution(A,b))
    {
        return 0;
    }
    return 1;
}

// Forward elimination
int forwardElimination(matrix *A, matrix *b)
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
int backwardSubstitution(matrix *A, matrix* b)
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
int assembleLoadVector(struct meshInfo meshInfoDb,struct boundary boundaryDb[],matrix *loadVector)
{
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        int nodeId = boundaryDb[i].nodeId;
        loadVector->mat[nodeId-1][0] = boundaryDb[i].value;
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
int assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb,struct element elementDb[],matrix elemMatrix[],matrix *globalMatrix)
{
    for (int elemId = 0; elemId < meshInfoDb.elementNum; elemId++)
    {
        int leftNodePos = elementDb[elemId].nodeId[0]-1;
        int rightNodePos = elementDb[elemId].nodeId[1]-1;
        
        globalMatrix->mat[leftNodePos][leftNodePos] += elemMatrix[elemId].mat[0][0];
        globalMatrix->mat[rightNodePos][rightNodePos] += elemMatrix[elemId].mat[1][1];
        globalMatrix->mat[leftNodePos][rightNodePos] += elemMatrix[elemId].mat[0][1];
        globalMatrix->mat[rightNodePos][leftNodePos] += elemMatrix[elemId].mat[1][0];
    }
    return 1;
}
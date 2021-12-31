#include <stdio.h>
#include <malloc.h>
#include "1DFDM.h"
#include <string.h>
int main()
{
    // reading mesh file
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

    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[0].id),&(boundaryDb[0].value));
    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[1].id),&(boundaryDb[1].value));
    // end read mesh file

    // assemble element stiffness matrix
    




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
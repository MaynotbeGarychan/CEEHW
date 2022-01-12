#include <stdio.h>
#include <malloc.h>
#include "2DFEM.h"
#include <string.h>

int main()
{
    //      reading mesh file     //
    // -------------------------  //
    struct meshInfo meshInfoDb;
    FILE *fileIo = fopen("report2Mesh1.txt","rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    char readType[4];
    fscanf(fileIo,"%s %d %d %d %d",readType,&(meshInfoDb.id),&(meshInfoDb.nodeNum),&(meshInfoDb.elementNum),&(meshInfoDb.boundaryNum));

    struct node nodeDb[10];
    struct element elementDb[10];
    struct boundary boundaryDb[10];

    fscanf(fileIo,"%s %d %lf %lf",readType,&(nodeDb[0].id),&(nodeDb[0].x),&(nodeDb[0].y));
    fscanf(fileIo,"%s %d %lf %lf",readType,&(nodeDb[1].id),&(nodeDb[1].x),&(nodeDb[1].y));
    fscanf(fileIo,"%s %d %lf %lf",readType,&(nodeDb[2].id),&(nodeDb[2].x),&(nodeDb[2].y));
    fscanf(fileIo,"%s %d %lf %lf",readType,&(nodeDb[3].id),&(nodeDb[3].x),&(nodeDb[3].y));
    fscanf(fileIo,"%s %d %lf %lf",readType,&(nodeDb[4].id),&(nodeDb[4].x),&(nodeDb[4].y));
    fscanf(fileIo,"%s %d %lf %lf",readType,&(nodeDb[5].id),&(nodeDb[5].x),&(nodeDb[5].y));

    fscanf(fileIo,"%s %d %d %d %d",readType,&(elementDb[0].id),&(elementDb[0].nodeId[0]),&(elementDb[0].nodeId[1]),&(elementDb[0].nodeId[2]));
    fscanf(fileIo,"%s %d %d %d %d",readType,&(elementDb[1].id),&(elementDb[1].nodeId[0]),&(elementDb[1].nodeId[1]),&(elementDb[1].nodeId[2]));
    fscanf(fileIo,"%s %d %d %d %d",readType,&(elementDb[2].id),&(elementDb[2].nodeId[0]),&(elementDb[2].nodeId[1]),&(elementDb[2].nodeId[2]));
    fscanf(fileIo,"%s %d %d %d %d",readType,&(elementDb[3].id),&(elementDb[3].nodeId[0]),&(elementDb[3].nodeId[1]),&(elementDb[3].nodeId[2]));

    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[0].nodeId),&(boundaryDb[0].value));
    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[1].nodeId),&(boundaryDb[0].value));
    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[2].nodeId),&(boundaryDb[0].value));
    fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[3].nodeId),&(boundaryDb[0].value));

    // translate the information //
    // ------------------------- //
    int bounNodeIdVec[10];
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        bounNodeIdVec[i] = boundaryDb[i].nodeId;
    }
    int dof = meshInfoDb.nodeNum - meshInfoDb.boundaryNum;
    int unknownNodeVec[10];
    int unknownNodeVecLen = 0;
    for (int i = 0; i < meshInfoDb.nodeNum; i++)
    {
        if (search(nodeDb[i].id,bounNodeIdVec,meshInfoDb.boundaryNum) == 0)
        {
            unknownNodeVec[unknownNodeVecLen] = nodeDb[i].id;
            unknownNodeVecLen++;
        }
    }
    
    // assemble related matrix //
    // ----------------------- //
    // assemble element matrix
    matrix elemMatrix[10];
    for(int i = 0; i < meshInfoDb.elementNum; i++)
    {
        allocateMatrix(&(elemMatrix[i]),3,3);
        assembleElementStiffnessMatrix(elementDb[i],nodeDb,&(elemMatrix[i]));
    }
    // assemble global stiffness matrix
    matrix globalMatrix;
    initilizeMatrix(&globalMatrix,meshInfoDb.nodeNum,meshInfoDb.nodeNum);
    assembleGlobalStiffnessMatrix(meshInfoDb,elementDb,elemMatrix,&globalMatrix);
    //printMatrix(&globalMatrix);
    // make a system
    matrix loadVector;
    initilizeMatrix(&loadVector,meshInfoDb.nodeNum,1);
    matrix systemMatrix;
    newCombineMatrixCol(globalMatrix,loadVector,&systemMatrix);
    printMatrix(&systemMatrix);
    // delect unused rows
    matrix unknownSystemMatrix;
    allocateMatrix(&unknownSystemMatrix,dof,systemMatrix.numCol);
    deleteBoundaryRows(meshInfoDb,boundaryDb,systemMatrix,&unknownSystemMatrix);
    // apply boundary conditions
    matrix appliedSystemMatrix;
    applyBoundaryCondtion(unknownSystemMatrix,unknownNodeVec,meshInfoDb,boundaryDb,&appliedSystemMatrix);
    printMatrix(&appliedSystemMatrix);

    // solve
    //matrix idvec;
    //allocateMatrix(&idvec,dof,1);
    //gaussianEliminationFEM(&appliedSystemMatrix,&idvec);

    return 1;
}


/**
 * @brief matrix related function
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

// inverse the matrix
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

// get block of matrix
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

// round diagonal component
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

// forward elimination
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

// backward subtitution
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

// combine two matrix into one by column
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

// get an identity matrix
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

// tranpose the matrix
void transposeMatrix(matrix inputMat, matrix *outMat)
{
    allocateMatrix(outMat,inputMat.numCol,inputMat.numRow);
    for(int i = 0; i<outMat->numRow;i++)
    {
        for (int j = 0; j < outMat->numCol; j++)
        {
            outMat->mat[i][j] = inputMat.mat[j][i];
        }
    }
}

// copy matrix
void copyMatrix(matrix inputMat,matrix *outMat)
{
    allocateMatrix(outMat,inputMat.numRow,inputMat.numCol);
    for (int i = 0; i < outMat->numRow; i++)
    {
        for (int j = 0; j < outMat->numCol; j++)
        {
            outMat->mat[i][j] = inputMat.mat[i][j];
        }
    }
}

// calculate the det of 22 matrix
double calculateDetMatrix22(matrix T)
{
    double val = T.mat[0][0]*T.mat[1][1]-T.mat[0][1]*T.mat[1][0];
    return val;
}

// Direct method
int gaussianEliminationFEM(matrix *A, matrixInt *indexVec)
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

/**
 * @brief fem calculation related funtion
 * 
 */

void assembleElementStiffnessMatrix(struct element elementDb,struct node nodeDb[],matrix *elemMatrix)
{
    // calculate the transformation matrix
    matrix transMatrix;
    allocateMatrix(&transMatrix,2,2);
    assembleTransMatrix(elementDb,nodeDb,&transMatrix);
    // its related value
    double detJ = calculateDetMatrix22(transMatrix);
    printMatrix(&transMatrix);
    inverseMatrix(&transMatrix);
    matrix transInvtransMat;
    transposeMatrix(transMatrix,&transInvtransMat);

    // element stiffness matrix component
    matrix dfaidr;
    allocateMatrix(&dfaidr,3,2);
    dfaidr.mat[0][0] = -1;
    dfaidr.mat[0][1] = -1;
    dfaidr.mat[1][0] = 1;
    dfaidr.mat[1][1] = 0;
    dfaidr.mat[2][0] = 0;
    dfaidr.mat[2][1] = 1;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            elemMatrix->mat[i][j] = 0.5*detJ*funcElemMatrix(transInvtransMat,dfaidr,i,j);
        }
    }
    printMatrix(elemMatrix);
}

double funcElemMatrix(matrix tranInvJ, matrix dfaidr, int m, int n)
{
    double dr1dx = tranInvJ.mat[0][0];
    double dr2dx = tranInvJ.mat[0][1];
    double dr1dy = tranInvJ.mat[1][0];
    double dr2dy = tranInvJ.mat[1][1];

    double retVal = (dr1dx*dr1dx+dr1dy*dr1dy)*dfaidr.mat[m][0]*dfaidr.mat[n][0] + 
        (dr1dx*dr2dx+dr1dy*dr2dy)*dfaidr.mat[m][0]*dfaidr.mat[n][1] +
        (dr1dx*dr1dx+dr1dy*dr1dy)*dfaidr.mat[m][1]*dfaidr.mat[n][0] + 
        (dr2dx*dr2dx+dr2dy*dr2dy)*dfaidr.mat[m][1]*dfaidr.mat[n][1];
    return retVal;
}

void assembleTransMatrix(struct element elementDb, struct node nodeDb[], matrix *transMatrix)
{
    int nodeOnePos = elementDb.nodeId[0]-1;
    int nodeTwoPos = elementDb.nodeId[1]-1;
    int nodeThreePos = elementDb.nodeId[2]-1;

    transMatrix->mat[0][0] = -nodeDb[nodeOnePos].x+nodeDb[nodeTwoPos].x;
    transMatrix->mat[0][1] = -nodeDb[nodeOnePos].x+nodeDb[nodeThreePos].x;
    transMatrix->mat[1][0] = -nodeDb[nodeOnePos].y+nodeDb[nodeTwoPos].y;
    transMatrix->mat[1][1] = -nodeDb[nodeOnePos].y+nodeDb[nodeThreePos].y;
}

void assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb,struct element elementDb[],matrix elemMat[],matrix *globalMat)
{
    for (int idElem = 0; idElem < meshInfoDb.elementNum; idElem++)
    {
        for (int i = 0; i < 3; i++)
        {
            int globalRowPos = elementDb[idElem].nodeId[i]-1;
            for (int j = 0; j < 3; j++)
            {
                int globalColPos = elementDb[idElem].nodeId[j]-1;
                globalMat->mat[globalRowPos][globalColPos] = elemMat[idElem].mat[i][j];
            }
        }
    }
}

int deleteBoundaryRows(struct meshInfo meshInfoDb, struct boundary boundaryDb[], matrix inputMat, matrix *outMat)
{
    // boundary node id vector
    int bounNodeVec[10];
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        bounNodeVec[i] = boundaryDb[i].nodeId;
    }
    int row = 0;
    for (int i = 0; i < inputMat.numRow; i++)
    {
        if (search(i+1,bounNodeVec,meshInfoDb.boundaryNum) == 0)
        {
            for (int j = 0; j < outMat->numCol; j++)
            {
                outMat->mat[row][j] = inputMat.mat[i][j];
            }
            row++;
            //printMatrix(&outMat);
        }
    }
    int dof = meshInfoDb.nodeNum - meshInfoDb.boundaryNum;
    if (dof != row)
    {
        return 0;
    }
    return 1;
}

void applyBoundaryCondtion(matrix inputMat,int unkownNodeIdVec[], struct meshInfo meshInfoDb,struct boundary boundaryDb[], matrix *outMat)
{
    int lvPos = inputMat.numCol-1;
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        double val = boundaryDb[i].value;
        int pos = boundaryDb[i].nodeId-1;
        for (int j = 0; j < inputMat.numRow; j++)
        {
            inputMat.mat[j][lvPos] -= inputMat.mat[j][pos]*val;
            inputMat.mat[j][pos] = 0.0;
        }
    }
    allocateMatrix(outMat,inputMat.numRow,inputMat.numRow+1);
    for (int i = 0; i < outMat->numRow; i++)
    {
        outMat->mat[i][outMat->numCol-1] = inputMat.mat[i][lvPos];
    }
    for (int i = 0; i < meshInfoDb.nodeNum-meshInfoDb.boundaryNum; i++)
    {
        for (int j = 0; j < outMat->numRow; j++)
        {
            int unknownNodePos = unkownNodeIdVec[i]-1;
            outMat->mat[j][i] = inputMat.mat[j][unknownNodePos];
        }
    }
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

int search(int val, int vec[], int vecLen)
{
    for (int i = 0; i < vecLen; i++)
    {
        if (val == vec[i])
        {
            return 1;
        }
    }
    return 0;
}


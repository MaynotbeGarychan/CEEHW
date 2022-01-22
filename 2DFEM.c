#include <stdio.h>
#include <malloc.h>
#include "2DFEM.h"
#include <string.h>
#include <math.h>

//#define OptionProblem
#define testCG

int main()
{
#ifndef testCG
    //      reading mesh file     //
    // -------------------------  //
    struct meshInfo meshInfoDb;
    // please specify the input file here
    FILE *fileIo = fopen("report2Mesh2.txt","rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    char readType[4];
    // reading overall information for mesh 
    fscanf(fileIo,"%s %d %d %d %d",readType,&(meshInfoDb.id),&(meshInfoDb.nodeNum),&(meshInfoDb.elementNum),&(meshInfoDb.boundaryNum));

    struct node nodeDb[1000];
    struct element elementDb[1000];
    struct boundary boundaryDb[1000];
    // read nodes information
    for (int i = 0; i < meshInfoDb.nodeNum; i++)
    {
        fscanf(fileIo,"%s %d %lf %lf",readType,&(nodeDb[i].id),&(nodeDb[i].x),&(nodeDb[i].y));
    }
    // read element information
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        fscanf(fileIo,"%s %d %d %d %d",readType,&(elementDb[i].id),&(elementDb[i].nodeId[0]),&(elementDb[i].nodeId[1]),&(elementDb[i].nodeId[2]));
    }
    // read boundary information
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        fscanf(fileIo,"%s %d %lf",readType,&(boundaryDb[i].nodeId),&(boundaryDb[i].value));
    }
    fclose(fileIo);

    // translate the information //
    // ------------------------- //
    int bounNodeIdVec[1000];
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        bounNodeIdVec[i] = boundaryDb[i].nodeId;
    }
    int dof = meshInfoDb.nodeNum - meshInfoDb.boundaryNum;
    int unknownNodeVec[1000];
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
    matrix elemMatrix[1000];
    for(int i = 0; i < meshInfoDb.elementNum; i++)
    {
        allocateMatrix(&(elemMatrix[i]),3,3);
        assembleElementStiffnessMatrix(elementDb[i],nodeDb,&(elemMatrix[i]));
    }
    // assemble global stiffness matrix
    matrix globalMatrix;
    initilizeMatrix(&globalMatrix,meshInfoDb.nodeNum,meshInfoDb.nodeNum);
    assembleGlobalStiffnessMatrix(meshInfoDb,elementDb,elemMatrix,&globalMatrix);
    printMatrix(&globalMatrix);
    // assemble loadvector
    matrix loadVector;
    initilizeMatrix(&loadVector,meshInfoDb.nodeNum,1);
#ifdef OptionProblem
    assembleLoadVector(meshInfoDb,elementDb,6,&loadVector);
#endif
    // make a system
    matrix systemMatrix;
    newCombineMatrixCol(globalMatrix,loadVector,&systemMatrix);
    printMatrix(&systemMatrix);

    // begin to solve the matrix //
    // ------------------------- //
    // delect unused rows
    matrixInt idvec;
    allocateMatrixInt(&idvec,meshInfoDb.nodeNum,1);
    for (int i = 0; i < idvec.numRow; i++)
    {
        idvec.mat[i][0] = nodeDb[i].id;
    }
    matrix unknownSystemMatrix;
    allocateMatrix(&unknownSystemMatrix,dof,systemMatrix.numCol);
    matrixInt unknownIdVec;
    allocateMatrixInt(&unknownIdVec,dof,1);
    deleteBoundaryRows(meshInfoDb,boundaryDb,systemMatrix,idvec,&unknownSystemMatrix,&unknownIdVec);
    // apply boundary conditions
    matrix appliedSystemMatrix;
    applyBoundaryCondtion(unknownSystemMatrix,unknownNodeVec,meshInfoDb,boundaryDb,&appliedSystemMatrix);
    printMatrix(&appliedSystemMatrix);
    // solve
    gaussianEliminationFEM(&appliedSystemMatrix,&unknownIdVec);
    printMatrix(&appliedSystemMatrix);

    // Output the result //
    // ----------------- //
    FILE *OutputIo = fopen("output.txt","w");
    if (OutputIo == NULL)
    {
        return 0;
    }
    for (int i = 0; i < appliedSystemMatrix.numRow; i++)
    {
        int nodeId = unknownIdVec.mat[i][0];
        fprintf(fileIo,"x%d %lf %lf %lf\n",nodeId,nodeDb[nodeId-1].x,nodeDb[nodeId-1].y,appliedSystemMatrix.mat[i][appliedSystemMatrix.numCol-1]);
    }
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        int nodeId = boundaryDb[i].nodeId;
        fprintf(fileIo,"x%d %lf %lf %lf\n",nodeId,nodeDb[nodeId-1].x,nodeDb[nodeId-1].y,boundaryDb[i].value);
    }
    fclose(OutputIo);
#endif
    matrix testMat;
    allocateMatrix(&testMat,5,6);
    createTestMatrixForCG(&testMat);
    printMatrix(&testMat);

    // gauss
    matrixInt idvec;
    allocateMatrixInt(&idvec,5,1);
    for (int i = 0; i < 5; i++)
    {
        idvec.mat[i][0] = i + 1;
    }
    gaussianEliminationFEM(&testMat,&idvec);
    //printMatrix(&testMat);
    for (int i = 0; i < testMat.numRow; i++)
    {
        printf("x%d: %lf\n",idvec.mat[i][0],testMat.mat[i][testMat.numCol-1]);
    }

    // iterative
    createTestMatrixForCG(&testMat);
    printMatrix(&testMat);
    matrix result;
    allocateMatrix(&result,5,1);
    printf("input matrix is\n");
    printMatrix(&testMat);
    conjugateSolveMatrix(testMat,1e-2,&result);


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
void printMatrix(const matrix *T)
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
    initializeIdentityMatrix(&inverseA,A->numRow);
    matrix C;
    if (!newCombineMatrixCol(*A,inverseA,&C))
    {
        return 0;
    }
    forwardElimination(&C);
    backwardSubtitution(&C);
    roundDiagonalComponent(&C);
    getBlockOfMatrix(C,0,C.numRow-1,A->numCol,C.numCol-1,A);
    return 1;
}

// get block of matrix
void getBlockOfMatrix(matrix A, int beginRowPos, int endRowPos, int beginColPos ,int endColPos, matrix *block)
/*
*   Target: get a block from a matrix
*   Input: A, begin, end     Output: block (no need to init)
*/
{
    int blockNumRow = endRowPos-beginRowPos+1;
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
void initializeIdentityMatrix(matrix *A, int numRow)
{
    allocateMatrix(A,numRow,numRow);
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
                if (abs(A->mat[j][i]) > abs(A->mat[targetPos][i]))
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

int addtoMatrix(matrix inMat, matrix *outMat)
{
    if (inMat.numRow != outMat->numRow || inMat.numCol != outMat->numCol)
    {
        printf("addtoMatrix: size of matrix is wrong");
        return 0;
    }
    for (int i = 0; i < outMat->numRow; i++)
    {
        for (int j = 0; j < outMat->numCol; j++)
        {
            outMat->mat[i][j] += inMat.mat[i][j];
        }
    }
    return 1;
}

int addMatrix(matrix inMat1, matrix inMat2, matrix *outMat)
{
    if (inMat1.numRow != inMat2.numRow || inMat1.numCol != inMat2.numCol)
    {
        printf("addMatrix: size of matrix is wrong");
        return 0;
    }
    for (int i = 0; i < outMat->numRow; i++)
    {
        for (int j = 0; j < outMat->numCol; j++)
        {
            outMat->mat[i][j] = inMat1.mat[i][j] + inMat2.mat[i][j];
        }
    }
    return 1;
}

int minustoMatrix(matrix inMat, matrix *outMat)
{
    if (inMat.numRow != outMat->numRow || inMat.numCol != outMat->numCol)
    {
        printf("minustoMatrix: size of matrix is wrong");
        return 0;
    }
    for (int i = 0; i < outMat->numRow; i++)
    {
        for (int j = 0; j < outMat->numCol; j++)
        {
            outMat->mat[i][j] -= inMat.mat[i][j];
        }
    }
    return 1;
}

int minusMatrix(matrix inMat1, matrix inMat2, matrix *outMat)
/*
*   Target: outMat = inMat1 - inMat2
*   Note: you have to allocate the memory for outMat firstly
*/
{
    if (inMat1.numRow != inMat2.numRow || inMat1.numCol != inMat2.numCol)
    {
        return 0;
    }
    for (int i = 0; i < outMat->numRow; i++)
    {
        for (int j = 0; j < outMat->numCol; j++)
        {
            outMat->mat[i][j] = inMat1.mat[i][j] - inMat2.mat[i][j];
        }
    }
    return 1;
}

int innerProduct(matrix inMat1, matrix inMat2, matrix *outMat)
/*
*   Target: Calculate the inner product of two matrix,
*           Their size should be accurate
*   Note: allocate for outMat before using this function.
*/
{
    if (inMat1.numCol != inMat2.numRow)
    {
        return 0;
    }
    for (int i = 0; i < outMat->numRow; i++)
    {
        for (int j = 0; j < outMat->numCol; j++)
        {
            for (int k = 0; k < inMat1.numCol; k++)
            {
                outMat->mat[i][j] += inMat1.mat[i][k]*inMat2.mat[k][j];
            }
        }
    }
    return 1;
}

double dotProductVec(matrix vec1,matrix vec2)
{
    double val = 0;
    // col * col
    if (vec1.numCol == 1 && vec2.numCol ==1)
    {
        if (vec1.numRow != vec2.numRow)
        {
            printf("dimension of dot Product vector is wrong!");
            return 0;
        }
        for (int i = 0; i < vec1.numRow; i++)
        {
            val += vec1.mat[i][0]*vec2.mat[i][0];
        }
        return val;
    }
    // row * row
    if (vec1.numRow == 1 && vec2.numRow ==1)
    {
        if (vec1.numCol != vec2.numCol)
        {
            printf("dimension of dot Product vector is wrong!");
            return 0;
        }
        for (int i = 0; i < vec1.numCol; i++)
        {
            val += vec1.mat[0][i]*vec2.mat[0][i];
        }
        return val;
    }
    // col * row
    if (vec1.numCol == 1 && vec2.numRow == 1)
    {
        if (vec1.numRow != vec2.numCol)
        {
            printf("dimension of dot Product vector is wrong!");
            return 0;
        }
        for (int i = 0; i < vec1.numRow; i++)
        {
            val += vec1.mat[i][0]*vec2.mat[0][i];
        }
        return val;
    }
    // row * col
    if (vec1.numRow == 1 && vec2.numCol == 1)
    {
        if (vec1.numCol != vec2.numRow)
        {
            printf("dimension of dot Product vector is wrong!");
            return 0;
        }
        for (int i = 0; i < vec1.numCol; i++)
        {
            val += vec1.mat[0][i]*vec2.mat[i][0];
        }
        return val;
    }
    printf("dimension of dot Product vector is wrong!");
    return 0;
}

double normVector(matrix vec)
{
    double val = 0;
    if (vec.numRow == 1)
    {
        for (int i = 0; i < vec.numCol; i++)
        {
            val += vec.mat[0][i]*vec.mat[0][i];
        }
    }
    else if (vec.numCol == 1)
    {
        for (int i = 0; i < vec.numRow; i++)
        {
            val += vec.mat[i][0]*vec.mat[i][0];
        }
    }
    return sqrt(val);
}

void scaleMatrix(matrix inMat, double factor, matrix *outMat)
/*
*   Target: factor*inMat = outMat, 
*            Scale all the component of inMat,
*            and make a new outMat
*    Note: you need to init outMat before scaling
*/
{
    for (int i = 0; i < outMat->numRow; i++)
    {
        for (int j = 0; j < outMat->numCol; j++)
        {
            outMat->mat[i][j] = inMat.mat[i][j]*factor;
        }
    }
}

int copyMatrix(matrix inMat, matrix* outMat)
{
    if (inMat.numCol != outMat->numCol || inMat.numRow != outMat->numRow)
    {
        printf("void copyMatrix: size of matrix is wrong!");
        return 0;
    }
    for (int i = 0; i < inMat.numRow; i++)
    {
        for (int j = 0; j < inMat.numCol; j++)
        {
            outMat->mat[i][j] = inMat.mat[i][j];
        }
    }
    return 1;
}

void conjugateSolveMatrix(const matrix systemMatrix, double tolerance, matrix *result)
/*
*   Target: Conjuate iteration to solve matrix
*   Note: init the result vector before using it.
*/
{
    //                Initilization              //
    // ========================================  //
    // prepare the matrix and vector
    matrix A;
    matrix b;
    getBlockOfMatrix(systemMatrix,0,systemMatrix.numRow-1,
        0,systemMatrix.numCol-2,&A);
    getBlockOfMatrix(systemMatrix,0,systemMatrix.numRow-1,systemMatrix.numCol-1,
        systemMatrix.numCol-1,&b);
    int dof = systemMatrix.numRow;

    // initial guess: x <= x0
    matrix x, x0;
    allocateMatrix(&x,dof,1);
    allocateMatrix(&x0,dof,1);
    x0.mat[0][0] = 1.4;
    x0.mat[1][0] = -9;
    x0.mat[2][0] = 23;
    x0.mat[3][0] = 4;
    x0.mat[4][0] = 1.9;
    copyMatrix(x0,&x);
    freeMatrix(&x0);

    // calculation the initial residual
    matrix Ax;
    allocateMatrix(&Ax,dof,1);
    innerProduct(A,x,&Ax);
    matrix residual;
    allocateMatrix(&residual,dof,1);
    minusMatrix(b,Ax,&residual);
    freeMatrix(&Ax);

    // copy the residual to the search direction(p)
    matrix p;
    allocateMatrix(&p,dof,1);
    copyMatrix(residual,&p);

    //                Iteration                  //
    // ========================================  //
    double alpha, beta, ratio;
    matrix Ap;
    allocateMatrix(&Ap,dof,1);

    matrix temp, oldResidual;
    allocateMatrix(&temp,dof,1);
    allocateMatrix(&oldResidual,dof,1);

    for (int i = 0; i < 100; i++)
    {
        // calculate alpha
        innerProduct(A,p,&Ap);
        alpha = dotProductVec(residual,residual)/dotProductVec(p,Ap);
        //alpha = dotProductVec(p,residual)/dotProductVec(p,Ap);
        // calculate new x
        scaleMatrix(p,alpha,&temp);
        addtoMatrix(temp,&x);
        // calculate new residual
        scaleMatrix(Ap,alpha,&temp);
        copyMatrix(residual,&oldResidual);
        minustoMatrix(temp,&residual);

        // check convergence
        ratio = normVector(residual)/normVector(b);
        if (ratio < tolerance)
        {
            break;
        }
        
        // calculate beta
        beta = dotProductVec(residual,residual)/dotProductVec(oldResidual,oldResidual);
        //beta = dotProductVec(residual,Ap)/dotProductVec(p,Ap);
        scaleMatrix(p,beta,&temp);
        addMatrix(residual,temp,&p);
    }
    
    //                   Output                  //
    // ========================================  //
    /*
    if (normVector(residual) > tolerance)
    {
        printf("conjugate method exceed the limit of iteration!");
    }
    // give the x to the result
    for (int i = 0; i < x.numRow; i++)
    {
        result->mat[i][0] = x.mat[i][0];
    }
    printf("the result after CG iteration is:\n");
    printMatrix(result);
    */
}

// this is a function platform to create a matrix to test the CG
void createTestMatrixForCG(matrix* testMat)
{
    // diagonal 
    for (int i = 0; i < testMat->numRow; i++)
    {
        testMat->mat[i][i] = i+1;
    }
    // non diagonal part
    for (int i = 0; i < testMat->numRow; i++)
    {
        for (int j = i+1; j < testMat->numRow; j++)
        {
            testMat->mat[i][j] = i+1;
            testMat->mat[j][i] = i+1;
        }
    }
    // bvec
    testMat->mat[0][testMat->numCol-1] = 1;
    testMat->mat[1][testMat->numCol-1] = 25;
    testMat->mat[2][testMat->numCol-1] = 31;
    testMat->mat[3][testMat->numCol-1] = -4;
    testMat->mat[4][testMat->numCol-1] = 5;
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
    //printMatrix(&transMatrix);
    // its related value
    double detJ = calculateDetMatrix22(transMatrix);
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

    double x1sady1s = tranInvJ.mat[0][0]*tranInvJ.mat[0][0]+tranInvJ.mat[0][1]*tranInvJ.mat[0][1];
    double x1x2ady1y2 = tranInvJ.mat[0][0]*tranInvJ.mat[0][1]+tranInvJ.mat[1][0]*tranInvJ.mat[1][1];
    double y1sady2s = tranInvJ.mat[1][0]*tranInvJ.mat[1][0]+tranInvJ.mat[1][1]*tranInvJ.mat[1][1];

    double retVal = x1sady1s*dfaidr.mat[m][0]*dfaidr.mat[n][0]+
        x1x2ady1y2*(dfaidr.mat[m][0]*dfaidr.mat[n][1]+dfaidr.mat[m][1]*dfaidr.mat[n][0])+
        y1sady2s*dfaidr.mat[m][1]*dfaidr.mat[n][1];

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
                globalMat->mat[globalRowPos][globalColPos] += elemMat[idElem].mat[i][j];
            }
        }
        //printMatrix(globalMat);
    }
}

int deleteBoundaryRows(struct meshInfo meshInfoDb, struct boundary boundaryDb[], matrix inputMat, matrixInt idVec,
                    matrix *outMat, matrixInt *unknownIdVec)
{
    // boundary node id vector
    int bounNodeVec[1000];
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
            unknownIdVec->mat[row][0] = idVec.mat[i][0];
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

void assembleLoadVector(struct meshInfo meshInfoDb, struct element elementDb[], double RHSvalue, matrix *loadVector)
{
    double val = -RHSvalue/6;
    initilizeMatrix(loadVector,meshInfoDb.nodeNum,1);
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        loadVector->mat[elementDb[i].nodeId[0]-1][0] += val;
        loadVector->mat[elementDb[i].nodeId[1]-1][0] += val;
        loadVector->mat[elementDb[i].nodeId[2]-1][0] += val;
    }
}


/**
 * @brief math related functions
 * 
 */
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


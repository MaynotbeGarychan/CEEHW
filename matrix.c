/*
*		Matrix operation, matrix solver for FEM analysis
*       Header file: matrix.h
*		Author: CHEN Jiawei, the University of Tokyo
*		Date:	2022/01/26
*/

#pragma once
#include "matrix.h"
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void allocateMatrix(matrix* T, int numRow, int numCol)
{
    T->mat = (double**)malloc(numRow * sizeof(double*));
    for (int i = 0; i < numRow; i++)
    {
        T->mat[i] = (double*)malloc(numCol * sizeof(double));
    }
    T->numCol = numCol;
    T->numRow = numRow;
}

void allocateMatrixInt(matrixInt* T, int numRow, int numCol)
{
    T->mat = (int**)malloc(numRow * sizeof(int*));
    for (int i = 0; i < numRow; i++)
    {
        T->mat[i] = (int*)malloc(numCol * sizeof(int));
    }
    T->numCol = numCol;
    T->numRow = numRow;
}

// initialize the value in matrix
void initilizeMatrix(matrix* T, int numRow, int numCol)
{
    allocateMatrix(T, numRow, numCol);
    for (int i = 0; i < numRow; i++)
    {
        for (int j = 0; j < numCol; j++)
        {
            T->mat[i][j] = 0.0;
        }
    }
}

// zeros the value in matrix
void ZeroMatrix(matrix* T)
{
	for (int i = 0; i < T->numRow; i++)
	{
		for (int j = 0; j < T->numCol; j++)
		{
			T->mat[i][j] = 0.0;
		}
	}
}

// free the matrix
void freeMatrix(matrix* T)
{
    for (int i = 0; i < T->numCol; i++)
    {
        free(T->mat[i]);
    }
    free(T->mat);
}

// print the matrix
void printMatrix(const matrix* T)
{
    for (int i = 0; i < T->numRow; i++)
    {
        for (int j = 0; j < T->numCol; j++)
        {
            printf("%lf,", T->mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// inverse the matrix
int inverseMatrix(matrix* A)
{
    matrix inverseA;
    initializeIdentityMatrix(&inverseA, A->numRow);
    matrix C;
    if (!newCombineMatrixCol(*A, inverseA, &C))
    {
        return 0;
    }
    forwardElimination(&C);
    backwardSubtitution(&C);
    roundDiagonalComponent(&C);
    getBlockOfMatrix(C, 0, C.numRow - 1, A->numCol, C.numCol - 1, A);
    return 1;
}

// get block of matrix
void getBlockOfMatrix(matrix A, int beginRowPos, int endRowPos, int beginColPos, int endColPos, matrix* block)
/*
*   Target: get a block from a matrix
*   Input: A, begin, end     Output: block (no need to init)
*/
{
    int blockNumRow = endRowPos - beginRowPos + 1;
    int blockNumCol = endColPos - beginColPos + 1;
    allocateMatrix(block, blockNumRow, blockNumCol);
    for (int i = 0; i < blockNumRow; i++)
    {
        for (int j = 0; j < blockNumCol; j++)
        {
            block->mat[i][j] = A.mat[beginRowPos + i][beginColPos + j];
        }
    }
}

// round diagonal component
int roundDiagonalComponent(matrix* A)
{
    if (A->numRow > A->numCol)
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
            A->mat[i][j] = A->mat[i][j] / A->mat[i][i];
        }
        A->mat[i][i] = 1.0;
    }
    return 1;
}

// forward elimination
int forwardElimination(matrix* A)
{
    if (A->numCol < A->numRow)
    {
        return 0;
    }
    for (int i = 0; i < A->numRow; i++)
    {
        for (int j = i + 1; j < A->numRow; j++)
        {
            if (A->mat[j][i] == 0)
            {
                continue;
            }
            double ratio = A->mat[j][i] / A->mat[i][i];
            A->mat[j][i] = 0;
            for (int k = i + 1; k < A->numCol; k++)
            {
                A->mat[j][k] = A->mat[j][k] - ratio * A->mat[i][k];
            }
        }
    }
    return 1;
}

// backward subtitution
int backwardSubtitution(matrix* A)
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
        for (int j = i - 1; j >= 0; j--)
        {
            if (A->mat[j][i] == 0)
            {
                continue;
            }
            double ratio = A->mat[j][i] / A->mat[i][i];
            A->mat[j][i] = 0.0;
            for (int k = A->numRow; k < A->numCol; k++)
            {
                A->mat[j][k] = A->mat[j][k] - ratio * A->mat[i][k];
            }
        }
    }
    return 1;
}

// combine two matrix into one by column
int newCombineMatrixCol(matrix A, matrix B, matrix* C)
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
    allocateMatrix(C, A.numRow, A.numCol + B.numCol);
    for (int i = 0; i < A.numRow; i++)
    {
        for (int j = 0; j < A.numCol; j++)
        {
            C->mat[i][j] = A.mat[i][j];
        }
        for (int k = 0; k < B.numCol; k++)
        {
            C->mat[i][A.numCol + k] = B.mat[i][k];
        }
    }
    return 1;
}

// get an identity matrix
void initializeIdentityMatrix(matrix* A, int numRow)
{
    allocateMatrix(A, numRow, numRow);
    for (int i = 0; i < A->numRow; i++)
    {
        for (int j = 0; j < A->numCol; j++)
        {
            if (i == j)
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
void transposeMatrix(matrix inputMat, matrix* outMat)
{
    allocateMatrix(outMat, inputMat.numCol, inputMat.numRow);
    for (int i = 0; i < outMat->numRow; i++)
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
    double val = T.mat[0][0] * T.mat[1][1] - T.mat[0][1] * T.mat[1][0];
    return val;
}

// Direct method
int gaussianEliminationSolveMatrix(matrix A, matrixInt* indexVec, matrix* result)
/*begin
    Target: solve ax = b, A = { a | b }
    Step:
        1. Forward elimination
        2. Backward substitution
    Flag:   1: True  0:False
end*/
{
    // forward elimination
    if (!forwardElimintationPivot(&A, indexVec))
    {
        return 0;
    }
    //printMatrix(A);
    // backward substitution
    if (!backwardSubtitution(&A))
    {
        return 0;
    }
    //printMatrix(A);
    // round the diagonal component
    if (!roundDiagonalComponent(&A))
    {
        return 0;
    }
    // get the result to the result vec
    int resultPos = A.numCol - 1;
    for (int i = 0; i < A.numRow; i++)
    {
        result->mat[i][0] = A.mat[i][resultPos];
    }
	printf("result by Gaussian Pivot Elimination is\n");
	printMatrix(result);
    return 1;
}

int forwardElimintationPivot(matrix* A, matrixInt* indexVec)
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
            for (int j = i + 1; j < A->numRow - 1; j++)
            {
                if (abs(A->mat[j][i]) > abs(A->mat[targetPos][i]))
                {
                    targetPos = j;
                }
            }
            swapRowMatrix(A, targetPos, i);
            swapRowMatrixInt(indexVec, targetPos, i);
        }
        // if non zero, we can eliminate it
        for (int j = i + 1; j < A->numRow; j++)
        {
            if (A->mat[j][i] == 0)
            {
                continue;
            }
            double ratio = A->mat[j][i] / A->mat[i][i];
            A->mat[j][i] = 0;
            for (int k = i + 1; k < A->numCol; k++)
            {
                A->mat[j][k] = A->mat[j][k] - ratio * A->mat[i][k];
            }
        }
    }
    return 1;
}

void swapRowMatrix(matrix* A, int rowOnePos, int rowTwoPos)
{
    double temp;
    for (int i = 0; i < A->numCol; i++)
    {
        temp = A->mat[rowOnePos][i];
        A->mat[rowOnePos][i] = A->mat[rowTwoPos][i];
        A->mat[rowTwoPos][i] = temp;
    }
}

void swapRowMatrixInt(matrixInt* A, int rowOnePos, int rowTwoPos)
{
    int temp;
    for (int i = 0; i < A->numCol; i++)
    {
        temp = A->mat[rowOnePos][i];
        A->mat[rowOnePos][i] = A->mat[rowTwoPos][i];
        A->mat[rowTwoPos][i] = temp;
    }
}

int addtoMatrix(matrix inMat, matrix* outMat)
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

int addMatrix(matrix inMat1, matrix inMat2, matrix* outMat)
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

int minustoMatrix(matrix inMat, matrix* outMat)
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

int minusMatrix(matrix inMat1, matrix inMat2, matrix* outMat)
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

int innerProduct(matrix inMat1, matrix inMat2, matrix* outMat)
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
                outMat->mat[i][j] += inMat1.mat[i][k] * inMat2.mat[k][j];
            }
        }
    }
    return 1;
}

double dotProductVec(matrix vec1, matrix vec2)
{
    double val = 0;
    // col * col
    if (vec1.numCol == 1 && vec2.numCol == 1)
    {
        if (vec1.numRow != vec2.numRow)
        {
            printf("dimension of dot Product vector is wrong!");
            return 0;
        }
        for (int i = 0; i < vec1.numRow; i++)
        {
            val += vec1.mat[i][0] * vec2.mat[i][0];
        }
        return val;
    }
    // row * row
    if (vec1.numRow == 1 && vec2.numRow == 1)
    {
        if (vec1.numCol != vec2.numCol)
        {
            printf("dimension of dot Product vector is wrong!");
            return 0;
        }
        for (int i = 0; i < vec1.numCol; i++)
        {
            val += vec1.mat[0][i] * vec2.mat[0][i];
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
            val += vec1.mat[i][0] * vec2.mat[0][i];
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
            val += vec1.mat[0][i] * vec2.mat[i][0];
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
            val += vec.mat[0][i] * vec.mat[0][i];
        }
    }
    else if (vec.numCol == 1)
    {
        for (int i = 0; i < vec.numRow; i++)
        {
            val += vec.mat[i][0] * vec.mat[i][0];
        }
    }
    return sqrt(val);
}

void scaleMatrix(matrix inMat, double factor, matrix* outMat)
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
            outMat->mat[i][j] = inMat.mat[i][j] * factor;
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

void conjugateGradientSolveMatrix(const matrix systemMatrix, double tolerance, matrix* result)
/*
*   Target: Conjuate iteration to solve matrix
*   Note: init the result vector before using it.
*/
{
    //                Initialization              //
    // ========================================  //
    // prepare the matrix and vector
    matrix A;
    matrix b;
    getBlockOfMatrix(systemMatrix, 0, systemMatrix.numRow - 1,
        0, systemMatrix.numCol - 2, &A);
    getBlockOfMatrix(systemMatrix, 0, systemMatrix.numRow - 1, systemMatrix.numCol - 1,
        systemMatrix.numCol - 1, &b);
    int dof = systemMatrix.numRow;

    // initial guess: x <= x0
    matrix x, x0;
    initilizeMatrix(&x, dof, 1);
    initilizeMatrix(&x0, dof, 1);
    copyMatrix(x0, &x);
    freeMatrix(&x0);

    // calculation the initial residual
    matrix Ax;
    initilizeMatrix(&Ax, dof, 1);
    innerProduct(A, x, &Ax);
    matrix residual;
    initilizeMatrix(&residual, dof, 1);
    printMatrix(&b);
    printMatrix(&Ax);
    minusMatrix(b, Ax, &residual);
    printMatrix(&residual);
    freeMatrix(&Ax);

    // copy the residual to the search direction(p)
    matrix p;
    initilizeMatrix(&p, dof, 1);
    copyMatrix(residual, &p);

    //                Iteration                  //
    // ========================================  //
    double alpha, beta, ratio;
    matrix newX, newResidual;
    initilizeMatrix(&newX, dof, 1);
    initilizeMatrix(&newResidual, dof, 1);

    matrix Ap;
    initilizeMatrix(&Ap, dof, 1);

    matrix temp, oldResidual;
    initilizeMatrix(&temp, dof, 1);
    initilizeMatrix(&oldResidual, dof, 1);

    int iteration = 0;
    int maxIterNum = 100;
    for (int i = 0; i < maxIterNum; i++)
    {
        // calculate alpha
        innerProduct(A, p, &Ap);
        alpha = dotProductVec(residual, residual) / dotProductVec(p, Ap);
        // calculate new x
        scaleMatrix(p, alpha, &temp);
        addMatrix(x, temp, &newX);
        ZeroMatrix(&temp);
        // calculate new residual
        innerProduct(A, newX, &temp);
        minusMatrix(b, temp, &newResidual);
        ZeroMatrix(&temp);

        // check convergence
        ratio = normVector(newResidual) / normVector(b);
        if (ratio < tolerance)
        {
            iteration++;
			copyMatrix(x, &newX);
			copyMatrix(residual, &newResidual);
            printf("CG is converged, Yes!\n");
            printf("CG: iteration num is %d\n", iteration);
            break;
        }

        // calculate beta
        beta = dotProductVec(newResidual, newResidual) / dotProductVec(residual, residual);
        scaleMatrix(p, beta, &temp);
        addMatrix(newResidual, temp, &p);
        ZeroMatrix(&temp);
        // go to new iteration
        copyMatrix(newX, &x);
        copyMatrix(newResidual, &residual);
        iteration++;
    }

    //                   Output                  //
    // ========================================  //
    if (iteration > maxIterNum)
    {
        printf("conjugate method exceed the limit of iteration!\n");
    }
    // give the x to the result
    for (int i = 0; i < x.numRow; i++)
    {
        result->mat[i][0] = x.mat[i][0];
    }
    printf("the result after CG iteration is:\n");
    printMatrix(result);
}

// this is a function platform to create a matrix to test the CG
void createTestMatrixForCG(matrix* testMat)
{
    // diagonal 
    for (int i = 0; i < testMat->numRow; i++)
    {
        testMat->mat[i][i] = 10;
    }
    // non diagonal part
    for (int i = 0; i < testMat->numRow; i++)
    {
        for (int j = i + 1; j < testMat->numRow; j++)
        {
            testMat->mat[i][j] = i + 1;
            testMat->mat[j][i] = i + 1;
        }
    }
    // bvec
    testMat->mat[0][testMat->numCol - 1] = 1;
    testMat->mat[1][testMat->numCol - 1] = 25;
    testMat->mat[2][testMat->numCol - 1] = 31;
    testMat->mat[3][testMat->numCol - 1] = -4;
    testMat->mat[4][testMat->numCol - 1] = 5;
}
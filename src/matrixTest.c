/*
*		Test platform for matrix operation, matrix solver
*       Header file: test.h
*		Author: CHEN Jiawei, the University of Tokyo
*		Date:	2022/01/26
*/
#pragma once
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "matrix.h"

int matrixTest()
{
    matrix testMat;
    allocateMatrix(&testMat, 5, 6);
    createTestMatrixForCG(&testMat);
	printf("input matrix is\n");
	printMatrix(&testMat);

    // gauss
    matrixInt idvec;
    allocateMatrixInt(&idvec, 5, 1);
    for (int i = 0; i < 5; i++)
    {
        idvec.mat[i][0] = i + 1;
    }
    matrix result;
    allocateMatrix(&result, testMat.numRow, 1);
    gaussianEliminationSolveMatrix(testMat, &idvec,&result);

    // iterative
    createTestMatrixForCG(&testMat);
    printf("input matrix is\n");
    printMatrix(&testMat);
    allocateMatrix(&result, 5, 1);
    conjugateGradientSolveMatrix(testMat, 1e-2, &result);

    return 1;
}


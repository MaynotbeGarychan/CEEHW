#pragma once

typedef struct
{
    double** mat;
    int numRow, numCol;
}matrix;

typedef struct
{
    int** mat;
    int numRow, numCol;
}matrixInt;

void allocateMatrix(matrix* T, int numRow, int numCol);
void allocateMatrixInt(matrixInt* T, int numRow, int numCol);
void initilizeMatrix(matrix* T, int numRow, int numCol);
void freeMatrix(matrix* T);
void printMatrix(const matrix* T);
double calculateDetMatrix22(matrix T);
int inverseMatrix(matrix* A);
void initializeIdentityMatrix(matrix* A, int numRow);
int newCombineMatrixCol(matrix A, matrix B, matrix* C);
int forwardElimination(matrix* A);
int backwardSubtitution(matrix* A);
int roundDiagonalComponent(matrix* A);
void getBlockOfMatrix(matrix A, int beginRowPos, int endRowPos, int beginColPos, int endColPos, matrix* block);
void transposeMatrix(matrix inputMat, matrix* outMat);
int gaussianEliminationSolveMatrix(matrix A, matrixInt* indexVec, matrix* result);
int forwardElimintationPivot(matrix* A, matrixInt* indexVec);
void swapRowMatrix(matrix* A, int rowOnePos, int rowTwoPos);
void swapRowMatrixInt(matrixInt* A, int rowOnePos, int rowTwoPos);
int addtoMatrix(matrix inMat, matrix* outMat);
int addMatrix(matrix inMat1, matrix inMat2, matrix* outMat);
int minustoMatrix(matrix inMat, matrix* outMat);
int minusMatrix(matrix inMat1, matrix inMat2, matrix* outMat);
double dotProductVec(matrix vec1, matrix vec2);
int innerProduct(matrix inMat1, matrix inMat2, matrix* outMat);
double normVector(matrix vec);
void conjugateGradientSolveMatrix(const matrix systemMatrix, double tolerance, matrix* result);
void scaleMatrix(matrix inMat, double factor, matrix* outMat);
void createTestMatrixForCG(matrix* testMat);
int copyMatrix(matrix inMat, matrix* outMat);
void ZeroMatrix(matrix* T);

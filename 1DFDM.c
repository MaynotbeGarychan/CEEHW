#include <stdio.h>
#include <malloc.h>
#include "1DFDM.h"
#include <string.h>
int main()
{
    matrix A;
    matrix b;    
    initilizeMatrix(&A,3,3);
    initilizeMatrix(&b,3,1);
    fillMatrix33Test(&A);
    fillMatrix31Test(&b);
    backwardSubstitution(&A,&b);
    for (int i = 0; i < b.numRow; i++)
    {
        double val = b.mat[i][0];
        printf("%f",val);
    }
    
    return 1;
}


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
    T->mat[1][0] = 0;
    T->mat[1][1] = 1;
    T->mat[1][2] = 2;
    T->mat[2][0] = 0;
    T->mat[2][1] = 0;
    T->mat[2][2] = 1;
}
void fillMatrix31Test(matrix *T)
{
    T->mat[0][0] = 12;
    T->mat[1][0] = 8;
    T->mat[2][0] = 3;
}

// Direct method
/*begin
    Target: solve Ax = b
    Step:
        1. Forward elimination A
        2. Combine {A,b}
        3. Backward substitution
    Flag:   1: True  0:False
end*/
int gaussianElimination(matrix *A,matrix *b)
{
    if (A->numRow != b->numRow)
    {
        return 0;
    }
    


    return 1;
}

// Forward elimination
void forwardElimination(matrix *T)
{
    for (int i = 0; i < T->numRow; i++)
    {
        for (int j = i+1; j < T->numRow; j++)
        {
            if (T->mat[j][i] == 0)
            {
                continue;
            }
            double ratio = T->mat[j][i]/T->mat[i][i];
            T->mat[j][i] = 0;
            for (int k = i+1; k < T->numCol; k++)
            {
                T->mat[j][k] = T->mat[j][k] - ratio * T->mat[i][k];
            }
        }
    }
}

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
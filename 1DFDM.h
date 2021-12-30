
typedef struct
{
    double **mat;
    int numRow, numCol;
}matrix;

void allocateMatrix(matrix *T, int numRow, int numCol);
void initilizeMatrix(matrix *T, int numRow, int numCol);
void freeMatrix(matrix *T);
void forwardElimination(matrix *T);
void fillMatrix33Test(matrix *T);
int backwardSubstitution(matrix *A, matrix* b);
int gaussianElimination(matrix *A,matrix *b);
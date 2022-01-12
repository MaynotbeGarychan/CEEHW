/**
 * @brief Matrix related
 * 
 */
typedef struct
{
    double **mat;
    int numRow, numCol;
}matrix;

typedef struct
{
    int **mat;
    int numRow, numCol;
}matrixInt;

void allocateMatrix(matrix *T, int numRow, int numCol);
void allocateMatrixInt(matrixInt *T, int numRow, int numCol);
void initilizeMatrix(matrix *T, int numRow, int numCol);
void freeMatrix(matrix *T);
void printMatrix(matrix *T);
double calculateDetMatrix22(matrix T);
int inverseMatrix(matrix *A);
void initializeIdentityMatrix(matrix *A);
int newCombineMatrixCol(matrix A, matrix B, matrix *C);
int forwardElimination(matrix *A);
int backwardSubtitution(matrix *A);
int roundDiagonalComponent(matrix *A);
void getBlockOfMatrix(matrix A, int beginRowPos, int endRowPos, int beginColPos ,int endColPos, matrix *block);
void transposeMatrix(matrix inputMat, matrix *outMat);
int gaussianEliminationFEM(matrix *A, matrixInt *indexVec);
int forwardElimintationPivot(matrix *A, matrixInt *indexVec);
void swapRowMatrix(matrix *A,int rowOnePos,int rowTwoPos);
void swapRowMatrixInt(matrixInt *A,int rowOnePos,int rowTwoPos);
void copyMatrix(matrix inputMat,matrix *outMat);

/**
 * @brief Database, solver and pre post processing
 * 
 */
struct node
{
    int id;
    double x;
    double y;
    //double z;
};

struct element
{
    int id;
    int nodeId[3];
};

struct boundary
{
    int nodeId;
    double value;
};

struct meshInfo
{
    int id;
    int nodeNum;
    int elementNum;
    int boundaryNum;
};

void assembleTransMatrix(struct element elementDb, struct node nodeDb[], matrix *transMatrix);
double funcElemMatrix(matrix tranInvJ, matrix dfaidr, int m, int n);
void assembleElementStiffnessMatrix(struct element elementDb,struct node nodeDb[],matrix *elemMatrix);
void assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb,struct element elementDb[],matrix elemMat[],matrix *globalMat);
int deleteBoundaryRows(struct meshInfo meshInfoDb, struct boundary boundaryDb[], matrix inputMat, matrix *outMat);
void applyBoundaryCondtion(matrix inputMat,int unkownNodeIdVec[], struct meshInfo meshInfoDb,struct boundary boundaryDb[], matrix *outMat);

/**
 * @brief math related functions
 * 
 */
double mathAbs(double a);
int search(int val, int vec[], int vecLen);
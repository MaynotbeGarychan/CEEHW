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
void fillMatrix33Test(matrix *T);
void printMatrix(matrix *T);
int gaussianEliminationFDM(matrix *A, matrixInt *indexVec);
void fillMatrix31Test(matrix *T);
int transposeMatrix(matrix *A, int rowPosOne, int rowPosTwo);
int getRowOfMatrix(matrix A, int rowPos, matrix *row);
int putRowOfMatrix(matrix row, int rowPos, matrix *A);
int inverseMatrix(matrix *A);
int backwardSubtitution(matrix *A);
int forwardElimination(matrix *A);
int roundDiagonalComponent(matrix *A);
void getBlockOfMatrix(matrix A, int beginRowPos, int endRowPos, int beginColPos ,int endColPos, matrix *block);
double dotProduct(matrix a, matrix b);
void initializeIdentityMatrix(matrix *A);
int newCombineMatrixCol(matrix A, matrix B, matrix *C);
int forwardElimintationPivot(matrix *A, matrixInt *indexVec);
void swapRowMatrix(matrix *A,int rowOnePos,int rowTwoPos);
void swapRowMatrixInt(matrixInt *A,int rowOnePos,int rowTwoPos);

/**
 * @brief math related functions
 * 
 */
typedef struct
{
    int *val;
    int length;
}arrayInt;

double mathAbs(double a);
//int allocateArrayInt(int *array, int length);

/**
 * @brief Database, solver and pre post processing
 * 
 */
struct node
{
    int id;
    double x;
    //double y;
    //double z;
};

struct element
{
    int id;
    int nodeId[2]; // 0: left 1: right
    //matrix stiffMatrix;
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

int readMeshHeader(FILE *fileIo, struct meshInfo *meshInfoDb);
int readMeshFile(const char* fileName, struct meshInfo *meshInfoDb, struct node *nodeDb);

int assembleElementStiffnessMatrix(struct element elementDb,struct node nodeDb[],matrix *elemMatrix);
int assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb,struct boundary boundaryDb[],struct element elementDb[],
                struct node nodeDb[], matrix elemMatrix[],matrix *globalMatrix);
int assembleLoadVector(struct meshInfo meshInfoDb, struct element elementDb[], struct node nodeDb[], matrix *loadVector);
int solveFEMSystem(matrix LHSMatrix, matrix RHSVector, struct meshInfo meshInfoDb, struct boundary boundaryDb[],
                struct node nodeDb[],matrix* result);
int findSameElementNode(struct meshInfo meshInfoDb, struct element elementDb[],int nodeId);
void reorderNodeList(struct meshInfo meshInfoDb, struct boundary boundaryDb[],
        struct node nodeDb[], matrixInt *reorderList);
void applyBoundaryCondtion(matrix *linearSystem, struct meshInfo meshInfoDb,struct boundary boundaryDb[]);
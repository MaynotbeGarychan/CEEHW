/**
 * @brief Matrix related
 * 
 */
typedef struct
{
    double **mat;
    int numRow, numCol;
}matrix;

void allocateMatrix(matrix *T, int numRow, int numCol);
void initilizeMatrix(matrix *T, int numRow, int numCol);
void freeMatrix(matrix *T);
int gaussForwardElimination(matrix *A, matrix *b);
void fillMatrix33Test(matrix *T);
int gaussBackwardSubstitution(matrix *A, matrix* b);
int gaussianElimination(matrix *A,matrix *b);
void fillMatrix31Test(matrix *T);
int transposeMatrix(matrix *A, int rowPosOne, int rowPosTwo);
int getRowOfMatrix(matrix A, int rowPos, matrix *row);
int putRowOfMatrix(matrix row, int rowPos, matrix *A);
int inverseMatrix(matrix *A);

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
    int nodeNum;
    int elementNum;
    int boundaryNum;
};

int readMeshHeader(const char* fileName, struct meshInfo *meshInfoDb);
int readMeshFile(const char* fileName, struct meshInfo *meshInfoDb, struct node *nodeDb);

int assembleElementStiffnessMatrix(struct element elementDb,struct node nodeDb[],matrix *elemMatrix);
int assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb,struct boundary boundaryDb[],struct element elementDb[],
                struct node nodeDb[], matrix elemMatrix[],matrix *globalMatrix);
int assembleLoadVector(struct meshInfo meshInfoDb, struct element elementDb[], struct node nodeDb[], matrix *loadVector);
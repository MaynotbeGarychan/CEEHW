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
int forwardElimination(matrix *A, matrix *b);
void fillMatrix33Test(matrix *T);
int backwardSubstitution(matrix *A, matrix* b);
int gaussianElimination(matrix *A,matrix *b);
void fillMatrix31Test(matrix *T);

/**
 * @brief Database and pre post processing
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
    int id;
    double value;
};

struct meshInfo
{
    int nodeNum;
    int elementNum;
    int boundaryNum;
};

typedef struct
{
    int elementId;
    matrix mat;
}elementMatrix;

typedef struct
{
    matrix mat;
    int dof;
}globalMatrix;

int readMeshHeader(const char* fileName, struct meshInfo *meshInfoDb);
int readMeshFile(const char* fileName, struct meshInfo *meshInfoDb, struct node *nodeDb);
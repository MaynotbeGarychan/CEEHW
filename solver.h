#pragma once
#include "matrix.h"
#include "mesh.h"
#include "macro.h"

enum problem
{
	WAVE_PRO = 1, POIS_PRO
};

enum matrixSolver
{
	GAUSSPIVOT_SOLVER = 1,CG_SOLVER
};

struct matrixSolverParam
{
	int matrixSolverType;
	double tolerance; // only used for cg
};

struct timeIntegrationParam
{
	double beta;
	double startTime;
	double endTime;
	double stepLength;
	int stepNum;
};

typedef struct
{
	int dimension;
	int solveProblem;
	int dof[MAX_NUM_TIMESTEP];
	struct matrixSolverParam solverParam;
	int appliedAbsorbingBoundary; // 0: no, 1: yes
	int boundaryNodeIdList[MAX_NUM_TIMESTEP][MAX_NUM_BOUD];
	int internalNodeIdList[MAX_NUM_TIMESTEP][MAX_NUM_NODE];

	int usedTimeInteScheme; // 0: no, 1: yes
	struct timeIntegrationParam timeInteParam;
}analysis;

//void assembleTransMatrix(struct element elementDb, struct node nodeDb[], matrix* transMatrix);
//double funcElemMatrix(matrix tranInvJ, matrix dfaidr, int m, int n);
//void assembleElementStiffnessMatrix(struct element elementDb, struct node nodeDb[], matrix* elemMatrix);
//void assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb, struct element elementDb[], matrix elemMat[], matrix* globalMat);
int deleteBoundaryRows(int totalNumNodes, int internalNodeNum, int internalNodeIdList[], matrix linearSystem, matrixInt idArray, matrix* slimLinearSystem, matrixInt* slimIdArray);
void applyBoundaryCondtion(matrix inputMat, int unkownNodeIdVec[], struct meshInfo meshInfoDb, struct boundary boundaryDb[], matrix* outMat);
//void assembleLoadVector(struct meshInfo meshInfoDb, struct element elementDb[], double RHSvalue, matrix* loadVector);
int search(int val, int vec[], int vecLen);
void initializeIdArray(mesh meshDb, matrixInt* idArray);
void applyBoundaryConditionAndDeleteCols(int staticBoundaryNum, int dynamicBoundaryNumStep, int internalNodeNum,
	struct boundary staticBoundaryDb[], struct boundaryDynamic dynamicBoundaryDbStep[], int internalNodeListStep[],
	matrix slimLinearSys, matrix* finalLinearSys);
int callMatrixSolver(int matrixSolverType, matrix finalLinearSys, matrixInt slimIdArray,
	matrix* resultArr);
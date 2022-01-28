#pragma once
#include "matrix.h"
#include "mesh.h"

void assembleTransMatrix(struct element elementDb, struct node nodeDb[], matrix* transMatrix);
double funcElemMatrix(matrix tranInvJ, matrix dfaidr, int m, int n);
void assembleElementStiffnessMatrix(struct element elementDb, struct node nodeDb[], matrix* elemMatrix);
void assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb, struct element elementDb[], matrix elemMat[], matrix* globalMat);
int deleteBoundaryRows(struct meshInfo meshInfoDb, struct boundary boundaryDb[], matrix inputMat, matrixInt idVec,
    matrix* outMat, matrixInt* unknownIdVec);
void applyBoundaryCondtion(matrix inputMat, int unkownNodeIdVec[], struct meshInfo meshInfoDb, struct boundary boundaryDb[], matrix* outMat);
void assembleLoadVector(struct meshInfo meshInfoDb, struct element elementDb[], double RHSvalue, matrix* loadVector);
int search(int val, int vec[], int vecLen);
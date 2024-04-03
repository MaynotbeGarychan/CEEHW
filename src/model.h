/************************************************************
FileName: model.h
Author: Chen
Date: 2022/02/02
Description: header file for model.c
***********************************************************/
#pragma once

#include "mesh.h"
#include "matrix.h"
#include "solver.h"

/*	Global
* 
*/
void assembleGlobalStiffnessMatrix(struct meshInfo meshInfoDb, struct element elementDb[], matrix elemMat[], matrix* globalMat);
void assembleElemSrcVecToLinearSys(matrix elemSrcVec, matrix* linearSystem);

/*	1D Wave static
*
*/
void assembleElementStiffnessMatrix1DWave(struct element elem, struct node nodeDb[], matrix* elemStiffMat);
void assembleElementSourceVector1DWave(mesh meshDb, double pVal, matrix* elemSrcVec);
int assemble1DWaveStatic(mesh meshDb, analysis analysisInfo, matrix* linearSystem);

/*	2D Poisson static
*
*/
int assemble2DPoissonStatic(mesh meshDb, analysis analysisInfo, matrix* linearSystem);
void assembleElementStiffnessMatrix2DPoisson(mesh meshDb, int index, matrix* elemStiffMat);
double funcElemMatrix2DPoissonTriElem(matrix tranInvJ, matrix dfaidr, int m, int n);
void assembleTransMatrix2DPoissonTriElem(struct element elementDb, struct node nodeDb[], matrix* transMatrix);
void assembleElemSrcVec2DPoisson(mesh meshDb, double pVal, matrix* loadVector);

/*	1D Wave dynamic
*
*/
void addVelocityMatrixToGlobalMatrix1DWaveDynamic(const double constant, matrix velocity, matrix* linearSys);
void assembleElementMassMatrix1DWaveDynamic(struct element elem, struct node nodeDb[], matrix* elemStiffMat);
void addMassMatrixToGlobalMatrix1DWaveDynamic(matrix acceleration, matrix massMatrix, matrix* linearSys);
void assembleGlobalStiffnessMatrix1DWaveDynamic(struct meshInfo meshInfoDb, struct element elementDb[], struct node nodeDb[], matrix* linearSystem);

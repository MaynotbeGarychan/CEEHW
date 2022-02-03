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

void assembleGlobalStiffnessMatrix(mesh meshDb, matrix elemMat[], matrix* globalMat);
void assembleElemSrcVecToLinearSys(matrix elemSrcVec, matrix* linearSystem);

void assembleElementStiffnessMatrix1DWave(mesh meshDb, int index, matrix* elemStiffMat);
void assembleElementSourceVector1DWave(mesh meshDb, double pVal, matrix* elemSrcVec);
int assemble1DWaveStatic(mesh meshDb, analysis analysisInfo, matrix* linearSystem);

int assemble2DPoissonStatic(mesh meshDb, analysis analysisInfo, matrix* linearSystem);
void assembleElementStiffnessMatrix2DPoisson(mesh meshDb, int index, matrix* elemStiffMat);
double funcElemMatrix2DPoissonTriElem(matrix tranInvJ, matrix dfaidr, int m, int n);
void assembleTransMatrix2DPoissonTriElem(struct element elementDb, struct node nodeDb[], matrix* transMatrix);
void assembleElemSrcVec2DPoisson(mesh meshDb, double pVal, matrix* loadVector);



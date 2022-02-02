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
void assembleElementStiffnessMatrix1DWave(mesh meshDb, int index, matrix* elemStiffMat);
void assembleElementSourceVector1DWave(mesh meshDb, matrix* elemSrcVec);
int assemble1DWaveStatic(mesh meshDb, analysis analysisInfo, matrix* linearSystem);
/************************************************************
FileName: femMain.h
Author: Chen
Date: 2022/02/02
Description: header file of the femMain.c
***********************************************************/
#pragma once
#include "Io.h"
#include "mesh.h"
#include "solver.h"
#include "model.h"
#include "matrix.h"
#include "macro.h"
#include "femMain.h"



int femMain(Io ioInfo);
int fem2dPoissonStatic(Io ioInfo, mesh meshDb, analysis analysisInfo, result* resultDb);
int fem1dWaveStatic(Io ioInfo, mesh meshDb, analysis analysisInfo, result* resultDb);
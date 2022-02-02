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

struct nodeScalarResult
{
	struct node nodeDb;
	double val;
	int timeStep;
};

typedef struct
{
	struct nodeScalarResult nodeScalarResultDb[MAX_NUM_NODE];
}result;


int femMain(Io ioInfo);
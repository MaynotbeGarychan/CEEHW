/************************************************************
FileName: translate.h
Author: Chen
Date: 2022/02/02
Description: header file of the translator.c
***********************************************************/
#pragma once
#include "mesh.h"
#include "translator.h"

struct nodeScalarResultStep
{
	struct node nodeDb[MAX_NUM_NODE];
	double val[MAX_NUM_NODE];
	double time;
};

typedef struct
{
	struct nodeScalarResultStep nodeScalarResultDb[MAX_NUM_TIMESTEP];
}result;

void translator(mesh* meshDb, analysis* analysisInfo);
void getBoundaryNodeListStep(int step, const mesh meshDb, int* boundaryNodeList[]);
void getInternalNodeListStep(int step, const mesh meshDb, const int* boundaryNodeList[],
	int* internalNodeList[]);
int getElemNodeNum(int elemType);

void saveScalarResultStep(int step, double time, mesh meshDb, analysis analysisInfo, matrix resultArr, matrixInt slimIdArray, result* resultDb);
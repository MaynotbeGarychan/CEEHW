/************************************************************
FileName: translate.h
Author: Chen
Date: 2022/02/02
Description: header file of the translator.c
***********************************************************/
#pragma once
#include "mesh.h"

void translator(mesh* meshDb, analysis* analysisInfo);
void getBoundaryNodeList(const mesh meshDb, int* boundaryNodeList[]);
void getInternalNodeList(const mesh meshDb, const int* boundaryNodeList[],
	int* internalNodeList[]);
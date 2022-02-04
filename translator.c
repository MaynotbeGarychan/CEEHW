/************************************************************
FileName: translate.c
Author: Chen
Date: 2022/02/02
Description: Translate the input information to the 
			  parameters for solver
***********************************************************/
#pragma once
#include "mesh.h"
#include "solver.h"
#include "translator.h"

void translator(mesh *meshDb, analysis *analysisInfo)
{
	/*	Calculate some basic param for time integration scheme
	*
	*/
	if (analysisInfo->usedTimeInteScheme == 1)
	{
		analysisInfo->timeInteParam.stepNum = (analysisInfo->timeInteParam.endTime - analysisInfo->timeInteParam.startTime)
			/ analysisInfo->timeInteParam.stepLength;
	}
	else
	{
		analysisInfo->usedTimeInteScheme = 0;
		analysisInfo->timeInteParam.stepLength = 0;
		analysisInfo->timeInteParam.startTime = 0;
		analysisInfo->timeInteParam.endTime = 0;
		analysisInfo->timeInteParam.beta = 0;
		analysisInfo->timeInteParam.stepNum = 0;
	}


	/*	1. Get the boundary node list, and then calculate
	*	 the internal node list.
	*	2. Calculate the boundary number, degree of freedom per step.
	* 
	*/
	switch (analysisInfo->usedTimeInteScheme)
	{
	case 1:
		for (int i = 0; i < analysisInfo->timeInteParam.stepNum; i++)
		{
			meshDb->boundaryInfoDb.totalBoundaryNumStep[i] =
				meshDb->boundaryInfoDb.staticBoundaryNum + meshDb->boundaryInfoDb.dynamicBoundaryNumStep[i];
			analysisInfo->dof[i] = meshDb->meshInfoDb.nodeNum - meshDb->boundaryInfoDb.totalBoundaryNumStep[i];
			getBoundaryNodeListStep(i + 1, *meshDb, analysisInfo->boundaryNodeIdList[i]);
			getInternalNodeListStep(i + 1, *meshDb, analysisInfo->boundaryNodeIdList[i], analysisInfo->internalNodeIdList[i]);
		}
		break;
	case 0:
		meshDb->boundaryInfoDb.dynamicBoundaryNumStep[0] = 0;
		meshDb->boundaryInfoDb.totalBoundaryNumStep[0] =
			meshDb->boundaryInfoDb.staticBoundaryNum + meshDb->boundaryInfoDb.dynamicBoundaryNumStep[0];

		analysisInfo->dof[0] = meshDb->meshInfoDb.nodeNum - meshDb->boundaryInfoDb.totalBoundaryNumStep[0];
		
		getBoundaryNodeListStep(1, *meshDb, analysisInfo->boundaryNodeIdList[0]);
		getInternalNodeListStep(1, *meshDb, analysisInfo->boundaryNodeIdList[0], analysisInfo->internalNodeIdList[0]);
		break;
	default:
		printf("translator.c : please check whether the timeInteParam is set?\n");
		break;
	}


	/*	Calculate the node number of one element
	*
	*/
	meshDb->meshInfoDb.elemNodeNum = getElemNodeNum(meshDb->meshInfoDb.elemType);


}

void getBoundaryNodeListStep(int step,const mesh meshDb, int *boundaryNodeList[])
{
	step = step - 1;
	for (int i = 0; i < meshDb.boundaryInfoDb.staticBoundaryNum; i++)
	{
		boundaryNodeList[i] = meshDb.staticBoundaryDb[i].nodeId;
	}
	for (int i = meshDb.boundaryInfoDb.staticBoundaryNum; i < meshDb.boundaryInfoDb.totalBoundaryNumStep[step]; i++)
	{
		boundaryNodeList[i] = meshDb.dynamicBoundaryDb[step][i].nodeId;
	}
}

void getInternalNodeListStep(int step, const mesh meshDb, const int *boundaryNodeList[],
	int *internalNodeList[])
{
	step = step - 1;
	int tmpVecPos = 0;
	for (int i = 0; i < meshDb.meshInfoDb.nodeNum; i++)
	{
		if (search(meshDb.nodeDb[i].id, boundaryNodeList, meshDb.boundaryInfoDb.totalBoundaryNumStep[step]) == 0)
		{
			internalNodeList[tmpVecPos] = meshDb.nodeDb[i].id;
			tmpVecPos++;
		}
	}
}

int getElemNodeNum(int elemType)
{
	switch (elemType)
	{
	case LINE_ELEM:
		return 2;
	case TRI_ELEM:
		return 3;
	default:
		printf("getElemNodeNum: unsupported element type!");
		break;
	}
}

void saveScalarResultStep(int step, double time, mesh meshDb, analysis analysisInfo, matrix resultArr, matrixInt slimIdArray, result* resultDb)
{
	int count = 0;
	resultDb->nodeScalarResultDb[step].time = time;
	for (int i = 0; i < resultArr.numRow; i++)
	{
		resultDb->nodeScalarResultDb[step].nodeDb[count].id = meshDb.nodeDb[slimIdArray.mat[i][0] - 1].id;
		resultDb->nodeScalarResultDb[step].nodeDb[count].x = meshDb.nodeDb[slimIdArray.mat[i][0] - 1].x;
		resultDb->nodeScalarResultDb[step].nodeDb[count].y = meshDb.nodeDb[slimIdArray.mat[i][0] - 1].y;

		resultDb->nodeScalarResultDb[step].val[count] = resultArr.mat[i][0];
		count++;
	}

	for (int i = 0; i < meshDb.boundaryInfoDb.staticBoundaryNum; i++)
	{
		resultDb->nodeScalarResultDb[step].nodeDb[count].id = meshDb.nodeDb[meshDb.staticBoundaryDb[i].nodeId - 1].id;
		resultDb->nodeScalarResultDb[step].nodeDb[count].x = meshDb.nodeDb[meshDb.staticBoundaryDb[i].nodeId - 1].x;
		resultDb->nodeScalarResultDb[step].nodeDb[count].y = meshDb.nodeDb[meshDb.staticBoundaryDb[i].nodeId - 1].y;

		resultDb->nodeScalarResultDb[step].val[count] = meshDb.staticBoundaryDb[i].value;
		count++;
	}

	for (int i = 0; i < meshDb.boundaryInfoDb.dynamicBoundaryNumStep[step]; i++)
	{
		resultDb->nodeScalarResultDb[step].nodeDb[count].id = meshDb.nodeDb[meshDb.dynamicBoundaryDb[step][i].nodeId - 1].id;
		resultDb->nodeScalarResultDb[step].nodeDb[count].x = meshDb.nodeDb[meshDb.dynamicBoundaryDb[step][i].nodeId - 1].x;
		resultDb->nodeScalarResultDb[step].nodeDb[count].y = meshDb.nodeDb[meshDb.dynamicBoundaryDb[step][i].nodeId - 1].y;

		resultDb->nodeScalarResultDb[step].val[count] = meshDb.dynamicBoundaryDb[step][i].value;
		count++;
	}
}
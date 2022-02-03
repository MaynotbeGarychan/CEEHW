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
	/*	Get the boundary node list, and then calculate
	*	 the internal node list.
	*  
	*/
	meshDb->boundaryInfoDb.totalBoundaryNum =
		meshDb->boundaryInfoDb.staticBoundaryNum + meshDb->boundaryInfoDb.dynamicBoundaryNum;
	getBoundaryNodeList(*meshDb, analysisInfo->boundaryNodeIdList);
	getInternalNodeList(*meshDb, analysisInfo->boundaryNodeIdList, analysisInfo->internalNodeIdList);

	/*	Calculate the node number of one element
	*
	*/
	meshDb->meshInfoDb.elemNodeNum = getElemNodeNum(meshDb->meshInfoDb.elemType);

	/*	Calculate the degree of freedom
	*	
	*/
	analysisInfo->dof = meshDb->meshInfoDb.nodeNum - meshDb->boundaryInfoDb.totalBoundaryNum;

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
		analysisInfo->timeInteParam.stepLength = 0;
		analysisInfo->timeInteParam.startTime = 0;
		analysisInfo->timeInteParam.endTime = 0;
		analysisInfo->timeInteParam.beta = 0;
		analysisInfo->timeInteParam.stepNum = 0;
	}

}

void getBoundaryNodeList(const mesh meshDb, int *boundaryNodeList[])
{
	for (int i = 0; i < meshDb.boundaryInfoDb.staticBoundaryNum; i++)
	{
		boundaryNodeList[i] = meshDb.staticBoundaryDb[i].nodeId;
	}
	for (int i = meshDb.boundaryInfoDb.staticBoundaryNum; i < meshDb.boundaryInfoDb.totalBoundaryNum; i++)
	{
		boundaryNodeList[i] = meshDb.dynamicBoundaryDb[i].nodeId;
	}
}

void getInternalNodeList(const mesh meshDb, const int *boundaryNodeList[],
	int *internalNodeList[])
{
	int tmpVecPos = 0;
	for (int i = 0; i < meshDb.meshInfoDb.nodeNum; i++)
	{
		if (search(meshDb.nodeDb[i].id, boundaryNodeList, meshDb.boundaryInfoDb.totalBoundaryNum) == 0)
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

void saveScalarResultStep(int step, mesh meshDb, matrix resultArr, matrixInt slimIdArray, result* resultDb)
{
	for (int i = 0; i < resultArr.numRow; i++)
	{
		resultDb->nodeScalarResultDb[]

		resultDb->nodeScalarResultDb[i].timeStep = 0;
		resultDb->nodeScalarResultDb[i].nodeDb = meshDb.nodeDb[slimIdArray.mat[i][0] - 1];
		resultDb->nodeScalarResultDb[i].val = resultArr.mat[i][0];
	}
}
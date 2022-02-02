/************************************************************
FileName: Io.c
Author: Chen
Date: 2022/02/02
Description: the source file for FEM Input and Output
***********************************************************/
#pragma once
#include "Io.h"
#include "mesh.h"
#include "solver.h"

int readTxt(Io ioInfo, mesh *meshDb, analysis *analysisInfo)
{
	const FILE* fileIo = fopen(ioInfo.inputDir, "rt");
	if (fileIo == NULL)
	{
		return 0;
	}

    /*	read the basic information for mesh
    *  Format: mesh id dimension elemType nodeNum elemNum BoundaryNum 
    */
    readMeshInfo(fileIo, &meshDb->meshInfoDb);

    /*	read the node information for mesh
    *  Format: node id xCoordinate yCoordinate
    */
    switch (meshDb->meshInfoDb.dimension)
    {
    case 1:
		for (int i = 0; i < meshDb->meshInfoDb.nodeNum; i++)
		{
			readNode1D(fileIo, &meshDb->nodeDb[i]);
		}
        break;
    case 2:
		for (int i = 0; i < meshDb->meshInfoDb.nodeNum; i++)
		{
			readNode2D(fileIo, &meshDb->nodeDb[i]);
		}
        break;
    default:
        printf("Unsupported dimension now!\n");
        break;
    }

    /*	read the element information for mesh
    *  Format: elem id nodeList[]
    */
	switch (meshDb->meshInfoDb.elemType)
	{
    case LINE_ELEM:
		for (int i = 0; i < meshDb->meshInfoDb.elementNum; i++)
		{
            readElemLine(fileIo, &meshDb->elementDb[i]);
		}
        break;
    case TRI_ELEM:
		for (int i = 0; i < meshDb->meshInfoDb.elementNum; i++)
		{
			readElemTri(fileIo, &meshDb->elementDb[i]);
		}
        break;
	default:
        printf("Unsupported element type now!\n");
		break;
	}

    /*	read the boundary information
    *  Format: boudhead id staticBoundaryNum dynamicBoundaryNum
    */
    readBoundaryInfo(fileIo, &meshDb->boundaryInfoDb);
    for (int i = 0; i < meshDb->boundaryInfoDb.staticBoundaryNum; i++)
    {
        readBoundary(fileIo, &meshDb->staticBoundaryDb[i]);
    }
    for (int i = 0; i < meshDb->boundaryInfoDb.dynamicBoundaryNum; i++)
    {
        readDynamicBoundary(fileIo, &meshDb->dynamicBoundaryDb[i]);
    }

    /*	read the analysis information
    *  Format: analysis problemName usedTimeIntegration
    */
	readAnalysisLine(fileIo, analysisInfo);

	/*	read the analysis solver information
	*  Format: analysis solverName
	*/
	readAnalysisLine(fileIo, analysisInfo);

	/*	read the analysis time integration information
	*  Format: analysis solverName
	*/
	readAnalysisLine(fileIo, analysisInfo);

	fclose(fileIo);
	return 1;
}

void readMeshInfo(const FILE* fileIo, struct meshInfo* meshInfoDb)
{
    char* readType[4];
    fscanf(fileIo, "%s %d %d %d %d %d %d\n", readType, &(meshInfoDb->id), &(meshInfoDb->dimension), &(meshInfoDb->elemType),
        &(meshInfoDb->nodeNum), &(meshInfoDb->elementNum), &(meshInfoDb->boundaryNum));
}

void readNode1D(const FILE* fileIo, struct node* node)
{
	char* readType[4];
	fscanf(fileIo, "%s %d %lf %lf\n", readType, &(node->id), &(node->x));
}

void readNode2D(const FILE* fileIo, struct node* node)
{
    char* readType[4];
    fscanf(fileIo, "%s %d %lf %lf\n", readType, &(node->id), &(node->x), &(node->y));
}

void readElem(const FILE* fileIo, struct element* elem)
{
    char* readType[4];
    fscanf(fileIo, "%s %d %d %d %d\n", readType, &(elem->id), &(elem->nodeId[0]), &(elem->nodeId[1]), &(elem->nodeId[2]));
}

void readElemLine(const FILE* fileIo, struct element* elem)
{
	char* readType[4];
	fscanf(fileIo, "%s %d %d %d %d\n", readType, &(elem->id), &(elem->nodeId[0]), &(elem->nodeId[1]));
}

void readElemTri(const FILE* fileIo, struct element* elem)
{
	char* readType[4];
	fscanf(fileIo, "%s %d %d %d %d\n", readType, &(elem->id), &(elem->nodeId[0]), &(elem->nodeId[1]), &(elem->nodeId[2]));
}

void readBoundary(const FILE* fileIo, struct boundary* boundary)
{
    char* readType[4];
    fscanf(fileIo, "%s %d %lf\n", readType, &(boundary->nodeId), &(boundary->value));
}

void readBoundaryInfo(const FILE* fileIo, struct boundaryInfo* boundaryInfoDb)
{
	char* readType[8];
    fscanf(fileIo, "%s %d %d %d\n", readType, &boundaryInfoDb->id, &boundaryInfoDb->staticBoundaryNum, &boundaryInfoDb->dynamicBoundaryNum);;
}

void readStaticBoundary(const FILE* fileIo, struct boundary* boundary)
{
	char* readType[8];
	fscanf(fileIo, "%s %d %lf\n", readType, &(boundary->nodeId), &(boundary->value));
}

void readDynamicBoundary(const FILE* fileIo, struct boundaryDynamic* boundary)
{
	char* readType[8];
	fscanf(fileIo, "%s %d %lf %lf\n", readType, &(boundary->nodeId), &boundary->time, &(boundary->value));
}

void readAnalysisLine(const FILE* fileIo, analysis* analysisInfo)
{
	char* readType[8];
	char* fillType[8];
	fscanf(fileIo, "%s %s ", readType, fillType);
	if (strcmp(fillType,"probName") == 0)
	{
		char* solveProbStr[4];
		fscanf(fileIo, "%s\n", solveProbStr);
		analysisInfo->solveProblem = mapProbStrToInt(solveProbStr);
	}
	else if (strcmp(fillType, "solvName") == 0)
	{
		char* solverTypeStr[4];
		fscanf(fileIo, "%s\n", solverTypeStr);
		analysisInfo->solverParam.matrixSolverType = mapSolverStrToInt(solverTypeStr);
	}
	else if (strcmp(fillType, "timeInte") == 0)
	{
		fscanf(fileIo, "%d %lf %lf %lf %lf\n", &analysisInfo->usedTimeInteScheme,
			&analysisInfo->timeInteParam.startTime, &analysisInfo->timeInteParam.endTime,
			&analysisInfo->timeInteParam.stepLength, &analysisInfo->timeInteParam.beta);
	}
}

int mapProbStrToInt(char* solveProbStr[4])
{
	if (strcmp(solveProbStr,"wave") == 0)
	{
        return WAVE_PRO;
	}
	else if (strcmp(solveProbStr, "pois") == 0)
	{
        return POIS_PRO;
	}
	else
	{
        printf("Unsupported problem now!");
	}
    return -1;
}

int mapSolverStrToInt(char* solveProbStr[4])
{
	if (strcmp(solveProbStr, "cgis") == 0)
	{
		return CG_SOLVER;
	}
	else if (strcmp(solveProbStr, "gpes") == 0)
	{
		return GAUSSPIVOT_SOLVER;
	}
	else
	{
		printf("Unsupported solver now!");
	}
	return -1;
}
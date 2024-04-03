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
#include "translator.h"

/*	Read Io
*
*/
int readTxt(Io ioInfo, mesh *meshDb, analysis *analysisInfo)
{
	const FILE* fileIo = fopen(ioInfo.inputDir, "rt");
	if (fileIo == NULL)
	{
		printf("radTxt: Cannot get input Io!\n");
		return 0;
	}

    /*	read the basic information for mesh
    *  Format: mesh id dimension elemTypeStr nodeNum elemNum BoundaryNum 
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
	int DynamicBoundarySetNum = 0;
    readBoundaryInfo(fileIo, &meshDb->boundaryInfoDb, &DynamicBoundarySetNum);
    for (int i = 0; i < meshDb->boundaryInfoDb.staticBoundaryNum; i++)
    {
        readBoundary(fileIo, &meshDb->staticBoundaryDb[i]);
    }
    for (int i = 0; i < DynamicBoundarySetNum; i++)
    {
		readDynamicBoundaryStep(fileIo, &meshDb->dynamicBoundaryDb[i], &(meshDb->boundaryInfoDb.dynamicBoundaryNumStep[i]));
    }

    /*	read the analysis information
    *  Format: analysis problemName usedTimeIntegration
    */
	/*	read the analysis solver information
	*  Format: analysis solverName
	*/
	/*	read the analysis time integration information
	*  Format: analysis timeInte usedTimeInteScheme start end stepLen beta
	*/
	int analysisLineNum = 0;
	readAnalysisHead(fileIo, &analysisLineNum);
	for (int i = 0; i < analysisLineNum; i++)
	{
		readAnalysisLine(fileIo, analysisInfo);
	}

	fclose(fileIo);
	return 1;
}

void readMeshInfo(const FILE* fileIo, struct meshInfo* meshInfoDb)
{
    char* readType[4];
	char* elemTypeStr[4];
    fscanf(fileIo, "%s %d %d %s %d %d %d\n", readType, &(meshInfoDb->id), &(meshInfoDb->dimension), elemTypeStr,
        &(meshInfoDb->nodeNum), &(meshInfoDb->elementNum), &(meshInfoDb->boundaryNum));
	meshInfoDb->elemType = mapElemTypeStrToInt(elemTypeStr);
}

int mapElemTypeStrToInt(char* solveProbStr[4])
{
	if (strcmp(solveProbStr, "line") == 0)
	{
		return LINE_ELEM;
	}
	else if (strcmp(solveProbStr, "tria") == 0)
	{
		return TRI_ELEM;
	}
	else
	{
		printf("Unsupported solver now!");
	}
	return -1;
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

void readBoundaryInfo(const FILE* fileIo, struct boundaryInfo* boundaryInfoDb, int *inteTimeStep)
{
	char* readType[8];
    fscanf(fileIo, "%s %d %d %d\n", readType, &boundaryInfoDb->id, &boundaryInfoDb->staticBoundaryNum, inteTimeStep);
}

void readStaticBoundary(const FILE* fileIo, struct boundary* boundary)
{
	char* readType[8];
	fscanf(fileIo, "%s %d %lf\n", readType, &(boundary->nodeId), &(boundary->value));
}

void readDynamicBoundaryStep(const FILE* fileIo, struct boundaryDynamic* boundary, int* numBoudStep)
{
	char* readType[8];
	fscanf(fileIo, "%s %lf %d\n", readType, &boundary->time, numBoudStep);
	for (int i = 0; i < *numBoudStep; i++)
	{
		fscanf(fileIo, "%s %d %lf\n", readType, &(boundary->nodeId), &(boundary->value));
	}
}

void readAnalysisHead(const FILE* fileIo, int *analysisNum)
{
	char* readType[8];
	fscanf(fileIo, "%s %d\n ", readType, analysisNum);
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
		analysisInfo->usedTimeInteScheme = 1;
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

/*	Write Io
*
*/
void writeTxt(Io ioInfo, mesh meshDb, analysis analysisInfo, result resultDb)
{
	FILE* OutputIo = fopen(ioInfo.outputDir, "w");
	if (OutputIo == NULL)
	{
		return 0;
	}
	
	switch (analysisInfo.usedTimeInteScheme)
	{
	case 1:
		for (int step = 0; step < analysisInfo.timeInteParam.stepNum; step++)
		{
			fprintf(OutputIo, "step %d time %lf\n", step + 1, resultDb.nodeScalarResultDb[step].time);
			for (int index = 0; index < meshDb.meshInfoDb.nodeNum; index++)
			{
				fprintf(OutputIo, "x%d %lf %lf %lf\n", resultDb.nodeScalarResultDb[step].nodeDb[index].id,
					resultDb.nodeScalarResultDb[step].nodeDb[index].x, resultDb.nodeScalarResultDb[step].nodeDb[index].y,
					resultDb.nodeScalarResultDb[step].val[index]);
			}
		}
		break;
	case 0:
		for (int index = 0; index < meshDb.meshInfoDb.nodeNum; index++)
		{
			fprintf(OutputIo, "x%d %lf %lf %lf\n", resultDb.nodeScalarResultDb[0].nodeDb[index].id,
				resultDb.nodeScalarResultDb[0].nodeDb[index].x, resultDb.nodeScalarResultDb[0].nodeDb[index].y,
				resultDb.nodeScalarResultDb[0].val[index]);
		}
		break;
	default:
		printf("writeTxt: plz specify the scheme!\n");
		break;
	}
}

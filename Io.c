/*
*		IO function, operation for pre, post processing
*       Header file: matrix.h
*		Author: CHEN Jiawei, the University of Tokyo
*		Date:	2022/01/26
*/
#pragma once
#include "Io.h"
#include "mesh.h"
#include "solver.h"

int readTxt(Io ioInfo, mesh *meshDb, analysis *analysisInfo)
{
	FILE* fileIo = fopen(ioInfo.inputDir, "rt");
	if (fileIo == NULL)
	{
		return 0;
	}

    // 
    meshDb->meshInfoDb.id = 0;

	fclose(fileIo);
	return 1;
}

int readMesh(const char *fileName,
    struct meshInfo *meshInfoDb, struct element *elementDb[] ,struct node *nodeDb[],struct boundary *boundaryDb[])
{
    FILE* fileIo = fopen(fileName, "rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    char readType[4];
    // reading overall information for mesh 
    fscanf(fileIo, "%s %d %d %d %d", readType, &(meshInfoDb->id), &(meshInfoDb->nodeNum), &(meshInfoDb->elementNum), &(meshInfoDb->boundaryNum));

    // read nodes information
    for (int i = 0; i < meshInfoDb->nodeNum; i++)
    {
        fscanf(fileIo, "%s %d %lf %lf", readType, &(nodeDb[i]->id), &(nodeDb[i]->x), &(nodeDb[i]->y));
    }
    // read element information
    for (int i = 0; i < meshInfoDb->elementNum; i++)
    {
        fscanf(fileIo, "%s %d %d %d %d", readType, &(elementDb[i]->id), &(elementDb[i]->nodeId[0]), &(elementDb[i]->nodeId[1]), &(elementDb[i]->nodeId[2]));
    }
    // read boundary information
    for (int i = 0; i < meshInfoDb->boundaryNum; i++)
    {
        fscanf(fileIo, "%s %d %lf", readType, &(boundaryDb[i]->nodeId), &(boundaryDb[i]->value));
    }
    fclose(fileIo);
	return 1;
}

void readMeshInfo(const FILE* fileIo, struct meshInfo* meshInfoDb)
{
    char* readType[4];
    fscanf(fileIo, "%s %d %d %d %d", readType,
        &(meshInfoDb->id), &(meshInfoDb->nodeNum), &(meshInfoDb->elementNum), &(meshInfoDb->boundaryNum));
}

void readNode(const FILE* fileIo, struct node* node)
{
    char* readType[4];
    fscanf(fileIo, "%s %d %lf %lf", readType, &(node->id), &(node->x), &(node->y));
}

void readElem(const FILE* fileIo, struct element* elem)
{
    char* readType[4];
    fscanf(fileIo, "%s %d %d %d %d", readType, &(elem->id), &(elem->nodeId[0]), &(elem->nodeId[1]), &(elem->nodeId[2]));
}

void readBoundary(const FILE* fileIo, struct boundary* boundary)
{
    char* readType[4];
    fscanf(fileIo, "%s %d %lf", readType, &(boundary->nodeId), &(boundary->value));
}
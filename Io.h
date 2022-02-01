#pragma once
#include "mesh.h"

typedef struct
{
	char* inputDir;
	char* outputDir;
}Io;

int readMesh(const char* fileName,
    struct meshInfo* meshInfoDb, struct element* elementDb[], struct node* nodeDb[], struct boundary* boundaryDb[]);


void readMeshInfo(const FILE* fileIo, struct meshInfo* meshInfoDb);
void readNode(const FILE* fileIo, struct node* node);
void readElem(const FILE* fileIo, struct element* elem);
void readBoundary(const FILE* fileIo, struct boundary* boundary);
#pragma once
#include "mesh.h"
#include "solver.h"
#include "translator.h"

typedef struct
{
	char* inputDir;
	char* outputDir;
}Io;

int readTxt(Io ioInfo, mesh* meshDb, analysis* analysisInfo);

void readMeshInfo(const FILE* fileIo, struct meshInfo* meshInfoDb);
void readNode1D(const FILE* fileIo, struct node* node);
void readNode2D(const FILE* fileIo, struct node* node);
void readElem(const FILE* fileIo, struct element* elem);
void readElemLine(const FILE* fileIo, struct element* elem);
void readElemTri(const FILE* fileIo, struct element* elem);
void readBoundary(const FILE* fileIo, struct boundary* boundary);
void readDynamicBoundaryStep(const FILE* fileIo, struct boundaryDynamic* boundary, int* numBoudStep);
void readStaticBoundary(const FILE* fileIo, struct boundary* boundary);
void readBoundaryInfo(const FILE* fileIo, struct boundaryInfo* boundaryInfoDb, int* inteTimeStep);
int mapSolverStrToInt(char* solveProbStr[4]);
int mapProbStrToInt(char* solveProbStr[4]);
void readAnalysisLine(const FILE* fileIo, analysis* analysisInfo);
int mapElemTypeStrToInt(char* solveProbStr[4]);
void readAnalysisHead(const FILE* fileIo, int* analysisNum);

void writeTxt(Io ioInfo, mesh meshDb, analysis analysisInfo, result resultDb);
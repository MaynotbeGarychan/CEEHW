/************************************************************
FileName: mesh.h
Author: Chen
Date: 2022/02/02
Description: header file of the mesh.c
***********************************************************/
#pragma once
#include <stdio.h>
#include "mesh.h"
#include "macro.h"

struct node
{
    int id;
    double x;
    double y;
    //double z;
};

enum elementType 
{
    LINE_ELEM = 1, TRI_ELEM
};

struct element
{
    int id;
    int nodeId[3];
};

struct meshInfo
{
    int id;
    int dimension;
    int elemType;
    int nodeNum;
    int elementNum;
    int boundaryNum;
};

struct boundary
{
	int nodeId;
	double value;
};

struct boundaryDynamic
{
    int nodeId;
    double time;
    double value;
};

struct boundaryInfo
{
    int id;
    int totalBoundaryNum;
    int staticBoundaryNum;
    int dynamicBoundaryNum;
};

typedef struct
{
	struct meshInfo meshInfoDb;
	struct element elementDb[MAX_NUM_ELEM];
	struct node nodeDb[MAX_NUM_NODE];
    struct boundaryInfo boundaryInfoDb;
    struct boundary staticBoundaryDb[MAX_NUM_BOUD];
    struct boundaryDynamic dynamicBoundaryDb[MAX_NUM_BOUD];
}mesh;

void addNode(int nodeId, double x, double y, struct node* nodeDb);
void addElem(int elemid, int nodeList[3], struct element* elemDb);
double calDistance(double x1, double y1, double x2, double y2);
int circleNode(int nodeId, struct element elem, struct node nodeDb[]);
int existInList(int index, int list[], int length);
int existinElemDb(int id, struct element elemDb[], int elemNum);
double calPolarAngle(double cartesianX, double cartesianY, double heartX, double heartY);
double maxNodeXCoor(struct node NodeDb[], int nodeNum);
double maxNodeYCoor(struct node NodeDb[], int nodeNum);
double minNodeXCoor(struct node NodeDb[], int nodeNum);
double minNodeYCoor(struct node NodeDb[], int nodeNum);

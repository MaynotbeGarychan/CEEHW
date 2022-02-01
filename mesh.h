#pragma once
#include <stdio.h>
#include "mesh.h"

#define MAX_NUM_ELEM 100
#define MAX_NUM_NODE 200
#define MAX_NUM_BOUD 100

struct node
{
    int id;
    double x;
    double y;
    //double z;
};

struct element
{
    int id;
    int nodeId[3];
};

struct boundary
{
    int nodeId;
    double value;
};

struct meshInfo
{
    int id;
    int nodeNum;
    int elementNum;
    int boundaryNum;
};

typedef struct
{
	struct meshInfo meshInfoDb;
	struct element elementDb[MAX_NUM_ELEM];
	struct node nodeDb[MAX_NUM_NODE];
    struct boundary boundaryDb[MAX_NUM_BOUD];
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

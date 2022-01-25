#pragma once
#include <stdio.h>

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

void addNode(int nodeId, double x, double y, struct node* nodeDb);
void addElem(int elemid, int nodeList[3], struct element* elemDb);
double calDistance(double x1, double y1, double x2, double y2);
int circleNode(int nodeId, struct element elem, struct node nodeDb[]);
int existInList(int index, int list[], int length);
int existinElemDb(int id, struct element elemDb[], int elemNum);
double calPolarAngle(double cartesianX, double cartesianY, double heartX, double heartY);
void orderPolygonNodeListByValRef(int* nodeList[], double valList[], int length);

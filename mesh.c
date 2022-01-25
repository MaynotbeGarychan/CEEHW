#pragma once
#include "mesh.h"
#include <math.h>

void addNode(int nodeId, double x, double y, struct node *nodeDb)
{
	nodeDb->id = nodeId;
	nodeDb->x = x;
	nodeDb->y = y;
}

void addElem(int elemid, int nodeList[3], struct element *elemDb)
{
	elemDb->id = elemid;
	elemDb->nodeId[0] = nodeList[0];
	elemDb->nodeId[1] = nodeList[1];
	elemDb->nodeId[2] = nodeList[2];
}

int circleNode(int nodeId, struct element elem, struct node nodeDb[])
/*
*	Target: check a node whether in a elem
*	Ret: 1: inside, 0: outside
*/
{
	// calculate heart
	double elemNode1X, elemNode2X, elemNode3X, elemNode1Y, elemNode2Y, elemNode3Y;
	elemNode1X = nodeDb[elem.nodeId[0] - 1].x;
	elemNode2X = nodeDb[elem.nodeId[1] - 1].x;
	elemNode3X = nodeDb[elem.nodeId[2] - 1].x;
	elemNode1Y = nodeDb[elem.nodeId[0] - 1].y;
	elemNode2Y = nodeDb[elem.nodeId[1] - 1].y;
	elemNode3Y = nodeDb[elem.nodeId[2] - 1].y;
	double heartX, heartY;
	heartX = (elemNode1X + elemNode2X + elemNode3X) / 3;
	heartY = (elemNode1Y + elemNode2Y + elemNode3Y) / 3;
	double addNodeX, addNodeY;
	addNodeX = nodeDb[nodeId - 1].x;
	addNodeY = nodeDb[nodeId - 1].y;

	// calculate diameter
	double diameter = 2 * calDistance(heartX, heartY, elemNode1X, elemNode1Y);

	// check if in the circle
	double val;
	for (int i = 0; i < 3; i++)
	{
		val = calDistance(addNodeX, addNodeY, nodeDb[elem.nodeId[i] - 1].x, nodeDb[elem.nodeId[i] - 1].y);
		if (val > diameter)
		{
			return 0;
		}
	}
	return 1;
}

double calDistance(double x1, double y1, double x2, double y2)
{
	double val = pow(x1 - x2, 2) + pow(y1 - y2, 2);
	return sqrt(val);
}

int existInList(int index, int list[], int length)
{
	for (int i = 0; i < length; i++)
	{
		if (list[i] == index)
		{
			return 1;
		}
	}
	return 0;
}

int existinElemDb(int id, struct element elemDb[], int elemNum)
{
	for (int i = 0; i < elemNum; i++)
	{
		if (elemDb[i].id == id)
		{
			return 1;
		}
	}
	return 0;
}

double calPolarAngle(double cartesianX, double cartesianY, double heartX, double heartY)
{
	double refX = cartesianX - heartX;
	double refY = cartesianY - heartY;

	double radius = calDistance(cartesianX, cartesianY, heartX, heartY);

	if (radius == 0)
	{
		printf("calPolarAngle : the point is just the heart, please enter other point!");
		return 0;
	}

	double angle;
	if (refY >= 0)
	{
		angle = acos(refX / radius);
	}
	else
	{
		angle = -acos(refX / radius);
	}
	return angle;
}

void orderPolygonNodeListByValRef(int *nodeList[], double valList[], int length)
{
	for (int i = 0; i < length; i++)
	{
		printf("%lf", nodeList[i]);
	}
}
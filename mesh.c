/*
*		Function, operation for Meshing, 2D Delauney triangulation
*       Header file: mesh.h
*		Author: CHEN Jiawei, the University of Tokyo
*		Date:	2022/01/26
*/
#pragma once
#include "mesh.h"
#define _USE_MATH_DEFINES
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
	double elemNode1X = nodeDb[elem.nodeId[0] - 1].x;
	double elemNode2X = nodeDb[elem.nodeId[1] - 1].x;
	double elemNode3X = nodeDb[elem.nodeId[2] - 1].x;
	double elemNode1Y = nodeDb[elem.nodeId[0] - 1].y;
	double elemNode2Y = nodeDb[elem.nodeId[1] - 1].y;
	double elemNode3Y = nodeDb[elem.nodeId[2] - 1].y;

	// give the val to the specific places, by the slope
	double x1, x2, x3, y1, y2, y3;
	if (elemNode1Y != elemNode2Y && elemNode1Y != elemNode3Y && elemNode2Y != elemNode3Y) // use anything is ok
	{
		x1 = elemNode1X;
		y1 = elemNode1Y;
		x2 = elemNode2X;
		y2 = elemNode2Y;
		x3 = elemNode3X;
		y3 = elemNode3Y;
	}
	else if (elemNode1Y == elemNode2Y) // use 13,23
	{
		x1 = elemNode1X;
		y1 = elemNode1Y;
		x2 = elemNode3X;
		y2 = elemNode3Y;
		x3 = elemNode2X;
		y3 = elemNode2Y;
	}
	else if (elemNode1Y == elemNode3Y) // use 12,23
	{
		x1 = elemNode1X;
		y1 = elemNode1Y;
		x2 = elemNode2X;
		y2 = elemNode2Y;
		x3 = elemNode3X;
		y3 = elemNode3Y;
	}
	else if (elemNode2Y == elemNode3Y) // use 12,13
	{
		x1 = elemNode2X;
		y1 = elemNode2Y;
		x2 = elemNode1X;
		y2 = elemNode1Y;
		x3 = elemNode3X;
		y3 = elemNode3Y;
	}
	else
	{
		printf("circleNode : the input circle may be wrong!");
	}

	double a1 = (y1 + y2) / 2;
	double a2 = (y2 + y3) / 2;
	double b1 = (x1 - x2) / (y2 - y1);
	double b2 = (x2 - x3) / (y3 - y2);
	double c1 = (x1 + x2) / 2;
	double c2 = (x2 + x3) / 2;

	double heartX = (b2 * c2 - b1 * c1 - a2 + a1) / (b2 - b1);
	double heartY = b1 * (heartX - c1) + a1;

	double addNodeX, addNodeY;
	addNodeX = nodeDb[nodeId - 1].x;
	addNodeY = nodeDb[nodeId - 1].y;

	// calculate radius
	double radius = calDistance(heartX, heartY, x1, y1);

	// check if in the circle
	double val = calDistance(addNodeX, addNodeY, heartX, heartY);
	if (val < radius)
	{
		return 1;
	}
	return 0;
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
		angle = 2 * M_PI - acos(refX / radius);
	}
	return angle;
}

double maxNodeXCoor(struct node NodeDb[], int nodeNum)
{
	int val = NodeDb[0].x;
	for (int i = 1; i < nodeNum; i++)
	{
		if (val < NodeDb[i].x)
		{
			val = NodeDb[i].x;
		}
	}
	return val;
}

double maxNodeYCoor(struct node NodeDb[], int nodeNum)
{
	int val = NodeDb[0].y;
	for (int i = 1; i < nodeNum; i++)
	{
		if (val < NodeDb[i].y)
		{
			val = NodeDb[i].y;
		}
	}
	return val;
}

double minNodeXCoor(struct node NodeDb[], int nodeNum)
{
	int val = NodeDb[0].x;
	for (int i = 1; i < nodeNum; i++)
	{
		if (val > NodeDb[i].x)
		{
			val = NodeDb[i].x;
		}
	}
	return val;
}

double minNodeYCoor(struct node NodeDb[], int nodeNum)
{
	int val = NodeDb[0].y;
	for (int i = 1; i < nodeNum; i++)
	{
		if (val > NodeDb[i].y)
		{
			val = NodeDb[i].y;
		}
	}
	return val;
}
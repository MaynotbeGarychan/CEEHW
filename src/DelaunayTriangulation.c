/*
*		2D Delauney triangulation to make mesh for FEM analysis
*       Header file: test.h
*		Author: CHEN Jiawei, the University of Tokyo 
*		Date:	2022/01/26
*/
#pragma once
#include "mesh.h"
#include "Io.h"

int meshDelauneyTest(Io ioInfo)
{
	// ----------------------------------------------------------  //
	//     reading seeds file, translate information               //
	// ----------------------------------------------------------  //

	struct node nodeDb[200] = {0};
	struct meshInfo meshInfoDb = { 0,0,0,0 };
	FILE* inputIo = fopen(ioInfo.inputDir, "rt");
	if (inputIo == NULL)
	{
		return 0;
	}
	readMeshInfo(inputIo, &meshInfoDb);
	for (int i = 0; i < meshInfoDb.nodeNum; i++)
	{
		readNode2D(inputIo, &nodeDb[i]);
	}
	fclose(inputIo);
	double maxNodeX = maxNodeXCoor(nodeDb, meshInfoDb.nodeNum);
	double minNodeX = minNodeXCoor(nodeDb, meshInfoDb.nodeNum);
	double maxNodeY = maxNodeYCoor(nodeDb, meshInfoDb.nodeNum);
	double minNodeY = minNodeYCoor(nodeDb, meshInfoDb.nodeNum);

	// ----------------------------------------------------------  //
	//     Begin to do Delaney Triangulation		               //
	// ----------------------------------------------------------  //

	// create an empty list to store triangle elements
	struct element elemDb[1000] = {0};

	// set rectangle which include all nodes
	//     -0.5 -0.5 1.5 1.5
	double xBound = (maxNodeX - minNodeX) * 0.2;
	double yBound = (maxNodeY - minNodeY) * 0.2;
	int boundNodeList[4] = { meshInfoDb.nodeNum + 1,meshInfoDb.nodeNum + 2,meshInfoDb.nodeNum + 3,meshInfoDb.nodeNum + 4 };
	addNode(boundNodeList[0], minNodeX - xBound, minNodeY - yBound, &nodeDb[boundNodeList[0] - 1]);
	addNode(boundNodeList[1], maxNodeX + xBound, minNodeY - yBound, &nodeDb[boundNodeList[1] - 1]);
	addNode(boundNodeList[2], maxNodeX + xBound, maxNodeY + yBound, &nodeDb[boundNodeList[2] - 1]);
	addNode(boundNodeList[3], minNodeX - xBound, maxNodeY + yBound, &nodeDb[boundNodeList[3] - 1]);

	// decompose rectangle by two triangle and add them to the triangle list
	meshInfoDb.elementNum = 2;
	int tempNodeList1[3] = { boundNodeList[0],boundNodeList[1],boundNodeList[3] };
	int tempNodeList2[3] = { boundNodeList[1],boundNodeList[2],boundNodeList[3] };
	addElem(1, tempNodeList1, &elemDb[0]);
	addElem(2, tempNodeList2, &elemDb[1]);

	// For each addutional node
	for (int nodeId = 1; nodeId <= meshInfoDb.nodeNum; nodeId++)
	{
		// check whether the node inside the element
		int elemCircleList[5] = {0};
 		int elemCircleNum = 0;
		for (int i = 0; i < meshInfoDb.elementNum; i++)
		{
			if (circleNode(nodeId, elemDb[i], &nodeDb))
			{
				elemCircleList[elemCircleNum] = elemDb[i].id;
				elemCircleNum++;
			}
		}

		// merge: delete, append all node list
		int polygonNodeList[100] = { 0 };
		int polygonNodeNum = 0;
		for (int i = 0; i < elemCircleNum; i++)
		{
			int deleElemId = elemCircleList[i];
			for (int j = 0; j < meshInfoDb.elementNum; j++)
			{
				// this is the dele Elem append their nodes
				if (elemDb[j].id == deleElemId)
				{
					for (int k = 0; k < 3; k++)
					{
						if (existInList(elemDb[j].nodeId[k],polygonNodeList,polygonNodeNum) == 0)
						{
							polygonNodeList[polygonNodeNum] = elemDb[j].nodeId[k];
							polygonNodeNum++;
						}
					}
					for (int g = j; g < meshInfoDb.elementNum; g++)
					{
						elemDb[g] = elemDb[g + 1];
					}
					meshInfoDb.elementNum--;
					break;
				}
			}
		}

		// make new element
		// calculate the angle of each polygon node
		double polygonNodeAngle[100] = { 0 };
		for (int i = 0; i < polygonNodeNum; i++)
		{
			polygonNodeAngle[i] = calPolarAngle(nodeDb[polygonNodeList[i] - 1].x, nodeDb[polygonNodeList[i] - 1].y,
				nodeDb[nodeId-1].x, nodeDb[nodeId-1].y);
		}

		// reorder the polygon node by a ref
		int reorderPolygonNodeList[100] = {0};
		for (int i = 0; i < polygonNodeNum; i++)
		{
			int id = polygonNodeList[i];
			double angle = polygonNodeAngle[i];
			int nodeRank = 1;
			for (int j = 0; j < polygonNodeNum; j++)
			{
				if (i != j)
				{
					double refAngle = polygonNodeAngle[j];
					if (angle > refAngle)
					{
						nodeRank++;
					}
				}
			}
			reorderPolygonNodeList[nodeRank-1] = id;
		}
		// also give the begin val to the end of reorder node list to make the final element
		reorderPolygonNodeList[polygonNodeNum] = reorderPolygonNodeList[0];
		// begin to make element
		for (int i = 0; i < polygonNodeNum; i++)
		{
			int nodeList[3] = { nodeId , reorderPolygonNodeList[i], reorderPolygonNodeList[i + 1] };
			for (int j = 1; j < 100; j++)
			{
				if (existinElemDb(j,elemDb,meshInfoDb.elementNum) == 0)
				{
					addElem(j, nodeList, &elemDb[meshInfoDb.elementNum]);
					meshInfoDb.elementNum++;
					break;
				}
			}
		}
	}

	// delete the bound element
	for (int i = 0; i < meshInfoDb.elementNum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (existInList(elemDb[i].nodeId[j], boundNodeList, 4))
			{
				for (int k = i; k < meshInfoDb.elementNum - 1; k++)
				{
					elemDb[k] = elemDb[k + 1];
				}
				meshInfoDb.elementNum--;
				i = 0; // if exise, after deket, restart the check
				break;
			}
		}
	}
	for (int j = 0; j < 3; j++)
	{
		if (existInList(elemDb[0].nodeId[j], boundNodeList, 4))
		{
			for (int k = 0; k < meshInfoDb.elementNum - 1; k++)
			{
				elemDb[k] = elemDb[k + 1];
			}
			meshInfoDb.elementNum--;
			break;
		}
	}


	// --------------------------------  //
	//     output the mesh to txt	     //
	// --------------------------------  //

	FILE* outputIo = fopen(ioInfo.outputDir, "w");
	if (outputIo == NULL)
	{
		return 0;
	}
	fprintf(outputIo, "mesh %d 1 2 tria %d %d 0\n", meshInfoDb.id, meshInfoDb.nodeNum, meshInfoDb.elementNum);
	for (int i = 0; i < meshInfoDb.nodeNum; i++)
	{
		fprintf(outputIo, "node %d %lf %lf\n", nodeDb[i].id, nodeDb[i].x, nodeDb[i].y);
	}
	for (int i = 0; i < meshInfoDb.elementNum; i++)
	{
		//fprintf(outputIo, "elem %d %d %d %d\n", elemDb[i].id, elemDb[i].nodeId[0], elemDb[i].nodeId[1], elemDb[i].nodeId[2]);
		fprintf(outputIo, "elem %d %d %d %d\n", i+1, elemDb[i].nodeId[0], elemDb[i].nodeId[1], elemDb[i].nodeId[2]);
	}
	fclose(outputIo);

	return 1;
}
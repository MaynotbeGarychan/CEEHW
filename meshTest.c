#pragma once
#include "mesh.h"
#include "Io.h"

int meshDelauneyTest(void)
{
	// set coordinates of nodes
	struct node nodeDb[100];
	FILE* inputIo = fopen("report3NodeSids.txt", "rt");
	if (inputIo == NULL)
	{
		return 0;
	}
	for (int i = 0; i < 6; i++)
	{
		readNode(inputIo, &nodeDb[i]);
	}
	fclose(inputIo);

	// create an empty list to store triangle elements
	struct element elemDb[100];

	// set rectangle which include all nodes
	//     -0.5 -0.5 1.5 1.5
	int boundNodeList[4] = { 7,8,9,10 };
	addNode(7, -0.5, -0.5, &nodeDb[6]);
	addNode(8, 1.5, -0.5, &nodeDb[7]);
	addNode(9, 1.5, 1.5, &nodeDb[8]);
	addNode(10, -0.5, 1.5, &nodeDb[9]);

	// decompose rectangle by two triangle and add them to the triangle list
	int elemNum = 2;
	int tempNodeList1[3] = { 7,8,10 };
	int tempNodeList2[3] = { 8,9,10 };
	addElem(1, tempNodeList1, &elemDb[0]);
	addElem(2, tempNodeList2, &elemDb[1]);

	// For each addutional node
	for (int nodeId = 1; nodeId <= 6; nodeId++)
	{
		// check whether the node inside the element
		int elemCircleList[5];
		int elemCircleNum = 0;
		for (int i = 0; i < elemNum; i++)
		{
			if (circleNode(nodeId, elemDb[i], &nodeDb))
			{
				elemCircleList[elemCircleNum] = elemDb[i].id;
				elemCircleNum++;
			}
		}

		// merge: delete, append all node list
		int polygonNodeList[20];
		int polygonNodeNum = 0;
		for (int i = 0; i < elemCircleNum; i++)
		{
			int deleElemId = elemCircleList[i];
			for (int j = 0; j < elemNum - 1; j++)
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
					for (int g = j; g < elemNum - 1; g++)
					{
						elemDb[g] = elemDb[g + 1];
					}
					elemNum--;
					break;
				}
			}
		}

		// make new element
		// reorder polygon node index
		// calculate the heart of the node list
		double heartX = 0;
		double heartY = 0;
		for (int i = 0; i < polygonNodeNum; i++)
		{
			heartX += nodeDb[polygonNodeList[i] - 1].x;
			heartY += nodeDb[polygonNodeList[i] - 1].y;
		}
		heartX = heartX / polygonNodeNum;
		heartY = heartY / polygonNodeNum;
		
		// calculate the angle of each polygon node
		double polygonNodeAngle[20];
		for (int i = 0; i < polygonNodeNum; i++)
		{
			polygonNodeAngle[i] = calPolarAngle(nodeDb[polygonNodeList[i] - 1].x, nodeDb[polygonNodeList[i] - 1].y,
				heartX, heartY);
		}

		// reorder the polygon node by a ref
		int reorderPolygonNodeList[20];
		for (int i = 0; i < polygonNodeNum; i++)
		{
			int id = polygonNodeList[i];
			double angle = polygonNodeAngle[i];
			int nodeRank = 1;
			for (int j = 0; j < polygonNodeNum; j++)
			{
				double refAngle = polygonNodeAngle[j];
				if (angle > refAngle)
				{
					nodeRank++;
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
				if (existinElemDb(j,elemDb,elemNum) == 0)
				{
					addElem(j, nodeList, &elemDb[elemNum]);
					elemNum++;
					break;
				}
			}
		}
	}

	// delete the bound element
	for (int i = 0; i < elemNum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (existInList(elemDb[i].nodeId[j], boundNodeList, 4))
			{
				for (int k = i; k < elemNum - 1; k++)
				{
					elemDb[k] = elemDb[k + 1];
				}
				elemNum--;
				i = 0;
				break;
			}
		}
	}
	for (int j = 0; j < 3; j++)
	{
		if (existInList(elemDb[0].nodeId[j], boundNodeList, 4))
		{
			for (int k = 0; k < elemNum - 1; k++)
			{
				elemDb[k] = elemDb[k + 1];
			}
			elemNum--;
			break;
		}
	}


	// output 
	FILE* outputIo = fopen("meshDealau.txt", "w");
	if (outputIo == NULL)
	{
		return 0;
	}
	for (int i = 0; i < 6; i++)
	{
		fprintf(outputIo, "node %d %lf %lf\n", nodeDb[i].id, nodeDb[i].x, nodeDb[i].y);
	}
	for (int i = 0; i < elemNum; i++)
	{
		fprintf(outputIo, "elem %d %d %d %d\n", elemDb[i].id, elemDb[i].nodeId[0], elemDb[i].nodeId[1], elemDb[i].nodeId[2]);
	}
	fclose(outputIo);

	return 1;
}
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
	//     -0.2 -0.2 1.2 1.2
	addNode(7, -0.2, -0.2, &nodeDb[6]);
	addNode(8, 1.2, -0.2, &nodeDb[7]);
	addNode(9, 1.2, 1.2, &nodeDb[8]);
	addNode(10, -0.2, 1.2, &nodeDb[9]);

	// decompose rectangle by two triangle and add them to the triangle list
	int elemNum = 4;
	int tempNodeList1[3] = { 1,7,8 };
	int tempNodeList2[3] = { 1,8,9 };
	int tempNodeList3[3] = { 1,9,10 };
	int tempNodeList4[3] = { 1,10,7 };
	addElem(1, tempNodeList1, &elemDb[0]);
	addElem(2, tempNodeList2, &elemDb[1]);
	addElem(3, tempNodeList3, &elemDb[2]);
	addElem(4, tempNodeList4, &elemDb[3]);

	// For each addutional node
	for (int nodeId = 2; nodeId <= 6; nodeId++)
	{
		// check whether the node inside the elemen
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

		// merge: delete, append all nodelist
		int polygonNodeList[20];
		int polygonNodeNum = 0;
		for (int i = 0; i < elemCircleNum; i++)
		{
			int deleElemId = elemCircleList[i];
			for (int j = 0; j < elemNum; j++)
			{
				// delete , re-arrange the node list
				if (elemDb[j].id > deleElemId)
				{
					elemDb[j-1] = elemDb[j];
				}
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
				}
			}
			elemNum--;
		}

		// make new element
		//reorder polygon node num

		// add extra one to the end
		for (int i = 0; i < polygonNodeNum; i++)
		{
			int nodeList[3] = { nodeId , polygonNodeList[i], polygonNodeList[i + 1] };
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
	
	
	// output 
	FILE* outputIo = fopen("meshDealau.txt", "w");
	if (outputIo == NULL)
	{
		return 0;
	}
	for (int i = 0; i < 2; i++)
	{
		fprintf(outputIo, "node %d %lf %lf\n", nodeDb[i].id, nodeDb[i].x, nodeDb[i].y);
	}
	for (int i = 0; i < 2; i++)
	{
		fprintf(outputIo, "elem %d %d %d %d\n", elemDb[i].id, elemDb[i].nodeId[0], elemDb[i].nodeId[1], elemDb[i].nodeId[2]);
	}
	fclose(outputIo);

	return 1;
}
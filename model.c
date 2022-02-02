/************************************************************
FileName: femMain.c
Author: Chen
Date: 2022/02/02
Description: src file for models
***********************************************************/
#pragma once

#include "mesh.h"
#include "matrix.h"
#include "solver.h"
#include "model.h"

/*  Global operation
*
*/
void assembleGlobalStiffnessMatrix(mesh meshDb, matrix elemMat[], matrix* globalMat)
{
	for (int idElem = 0; idElem < meshDb.meshInfoDb.elementNum; idElem++)
	{
		for (int i = 0; i < meshDb.meshInfoDb.elemNodeNum; i++)
		{
			int globalRowPos = meshDb.elementDb[idElem].nodeId[i] - 1;
			for (int j = 0; j < meshDb.meshInfoDb.elemNodeNum; j++)
			{
				int globalColPos = meshDb.elementDb[idElem].nodeId[j] - 1;
				globalMat->mat[globalRowPos][globalColPos] += elemMat[idElem].mat[i][j];
			}
		}
	}
	printf("Global stiffness matrix is: \n");
	printMatrix(globalMat);
}

/*	matrix and vector for 1D wave
*	G.E.: d2u/dx2 = 1
*/
int assemble1DWaveStatic(mesh meshDb, analysis analysisInfo, matrix* linearSystem)
{
	/*	Assemble element src vector
	*
	*/
	matrix elemSrcVec;
	initilizeMatrix(&elemSrcVec, meshDb.meshInfoDb.nodeNum, 1);
	assembleElementSourceVector1DWave(meshDb, &elemSrcVec);

	/*	Assemble element stiffness matrix
	*
	*/
	matrix elemStiffMat[MAX_NUM_ELEM];
	for (int i = 0; i < meshDb.meshInfoDb.elementNum; i++)
	{
		initilizeMatrix(&elemStiffMat[i], meshDb.meshInfoDb.elemNodeNum, meshDb.meshInfoDb.elemNodeNum);
		assembleElementStiffnessMatrix1DWave(meshDb, i, &elemStiffMat[i]);
	}

	/* Assemble global stiffness matrix
	*
	*/
	assembleGlobalStiffnessMatrix(meshDb, elemStiffMat, linearSystem);


	return 1;
}

void assembleElementStiffnessMatrix1DWave(mesh meshDb, int index, matrix* elemStiffMat)
/*	Note: initialize the matrix before enter this function
*			the index can be the list index 
*/
{
	// calculate the basic component
	double xleft = meshDb.nodeDb[meshDb.elementDb[index].nodeId[0] - 1].x;
	double xright = meshDb.nodeDb[meshDb.elementDb[index].nodeId[1] - 1].x;
	double k = 1 / (xright - xleft);
	// diagonal
	elemStiffMat->mat[0][0] = k;
	elemStiffMat->mat[1][1] = k;
	// non-diagonal
	elemStiffMat->mat[0][1] = -k;
	elemStiffMat->mat[1][0] = -k;
}

void assembleElementSourceVector1DWave(mesh meshDb, matrix* elemSrcVec)
{
	for (int i = 0; i < meshDb.meshInfoDb.elementNum; i++)
	{
		double temp = 0.5 * (meshDb.nodeDb[meshDb.elementDb[i].nodeId[0] - 1].x 
			- meshDb.nodeDb[meshDb.elementDb[i].nodeId[0] - 1].x);
		elemSrcVec->mat[meshDb.elementDb[i].nodeId[0] - 1][0] += temp;
		elemSrcVec->mat[meshDb.elementDb[i].nodeId[1] - 1][0] += temp;
	}
}

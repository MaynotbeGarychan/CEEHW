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

void assembleElemSrcVecToLinearSys(matrix elemSrcVec, matrix* linearSystem)
{
	for (int i = 0; i < linearSystem->numRow; i++)
	{
		linearSystem->mat[i][linearSystem->numCol - 1] += elemSrcVec.mat[i][0];
	}
}

/*	matrix and vector for 1D wave
*	G.E.: d2u/dx2 = p
*/
int assemble1DWaveStatic(mesh meshDb, analysis analysisInfo, matrix* linearSystem)
{
	/*	Assemble element src vector
	*
	*/
	matrix elemSrcVec;
	initilizeMatrix(&elemSrcVec, meshDb.meshInfoDb.nodeNum, 1);
	assembleElementSourceVector1DWave(meshDb, 1, &elemSrcVec);

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
	assembleElemSrcVecToLinearSys(elemSrcVec, linearSystem);
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

void assembleElementSourceVector1DWave(mesh meshDb, double pVal, matrix* elemSrcVec)
{
	for (int i = 0; i < meshDb.meshInfoDb.elementNum; i++)
	{
		double temp = 0.5 * (meshDb.nodeDb[meshDb.elementDb[i].nodeId[0] - 1].x 
			- meshDb.nodeDb[meshDb.elementDb[i].nodeId[1] - 1].x);
		elemSrcVec->mat[meshDb.elementDb[i].nodeId[0] - 1][0] += temp * pVal;
		elemSrcVec->mat[meshDb.elementDb[i].nodeId[1] - 1][0] += temp * pVal;
	}
}

/*	matrix and vector for 2D Poisson
*	G.E.: d2u/dx2 + d2u/dy2 = p
*/

int assemble2DPoissonStatic(mesh meshDb, analysis analysisInfo, matrix* linearSystem)
{
	/*	Assemble element stiffness matrix
	*
	*/
	matrix elemStiffMat[MAX_NUM_ELEM];
	for (int i = 0; i < meshDb.meshInfoDb.elementNum; i++)
	{
		initilizeMatrix(&elemStiffMat[i], meshDb.meshInfoDb.elemNodeNum, meshDb.meshInfoDb.elemNodeNum);
		assembleElementStiffnessMatrix2DPoisson(meshDb, i, &elemStiffMat[i]);
	}

	/*	Assemble element src vector
	*
	*/
	matrix elemSrcVec;
	initilizeMatrix(&elemSrcVec, meshDb.meshInfoDb.nodeNum, 1);
	assembleElemSrcVec2DPoisson(meshDb, 0, &elemSrcVec);

	/*	Assemble global stiffness matrix
	*
	*/
	assembleElemSrcVecToLinearSys(elemSrcVec, linearSystem);
	assembleGlobalStiffnessMatrix(meshDb, elemStiffMat, linearSystem);

	return 1;
}

void assembleElementStiffnessMatrix2DPoisson(mesh meshDb, int index, matrix* elemStiffMat)
{
	// calculate the transformation matrix
	matrix transMatrix;
	allocateMatrix(&transMatrix, 2, 2);
	assembleTransMatrix2DPoissonTriElem(meshDb.elementDb[index], meshDb.nodeDb, &transMatrix);
	// its related value
	double detJ = calculateDetMatrix22(transMatrix);
	inverseMatrix(&transMatrix);
	matrix transInvtransMat;
	transposeMatrix(transMatrix, &transInvtransMat);

	// element stiffness matrix component
	matrix dfaidr;
	allocateMatrix(&dfaidr, 3, 2);
	dfaidr.mat[0][0] = -1;
	dfaidr.mat[0][1] = -1;
	dfaidr.mat[1][0] = 1;
	dfaidr.mat[1][1] = 0;
	dfaidr.mat[2][0] = 0;
	dfaidr.mat[2][1] = 1;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			elemStiffMat->mat[i][j] = 0.5 * detJ * funcElemMatrix2DPoissonTriElem(transInvtransMat, dfaidr, i, j);
		}
	}
}

double funcElemMatrix2DPoissonTriElem(matrix tranInvJ, matrix dfaidr, int m, int n)
{
	double dr1dx = tranInvJ.mat[0][0];
	double dr2dx = tranInvJ.mat[0][1];
	double dr1dy = tranInvJ.mat[1][0];
	double dr2dy = tranInvJ.mat[1][1];

	double x1sady1s = tranInvJ.mat[0][0] * tranInvJ.mat[0][0] + tranInvJ.mat[0][1] * tranInvJ.mat[0][1];
	double x1x2ady1y2 = tranInvJ.mat[0][0] * tranInvJ.mat[0][1] + tranInvJ.mat[1][0] * tranInvJ.mat[1][1];
	double y1sady2s = tranInvJ.mat[1][0] * tranInvJ.mat[1][0] + tranInvJ.mat[1][1] * tranInvJ.mat[1][1];

	double retVal = x1sady1s * dfaidr.mat[m][0] * dfaidr.mat[n][0] +
		x1x2ady1y2 * (dfaidr.mat[m][0] * dfaidr.mat[n][1] + dfaidr.mat[m][1] * dfaidr.mat[n][0]) +
		y1sady2s * dfaidr.mat[m][1] * dfaidr.mat[n][1];

	return retVal;
}

void assembleTransMatrix2DPoissonTriElem(struct element elementDb, struct node nodeDb[], matrix* transMatrix)
{
	int nodeOnePos = elementDb.nodeId[0] - 1;
	int nodeTwoPos = elementDb.nodeId[1] - 1;
	int nodeThreePos = elementDb.nodeId[2] - 1;

	transMatrix->mat[0][0] = -nodeDb[nodeOnePos].x + nodeDb[nodeTwoPos].x;
	transMatrix->mat[0][1] = -nodeDb[nodeOnePos].x + nodeDb[nodeThreePos].x;
	transMatrix->mat[1][0] = -nodeDb[nodeOnePos].y + nodeDb[nodeTwoPos].y;
	transMatrix->mat[1][1] = -nodeDb[nodeOnePos].y + nodeDb[nodeThreePos].y;
}

void assembleElemSrcVec2DPoisson(mesh meshDb, double pVal, matrix* loadVector)
{
	double val = -pVal / 6;
	for (int i = 0; i < meshDb.meshInfoDb.elementNum; i++)
	{
		loadVector->mat[meshDb.elementDb[i].nodeId[0] - 1][0] += val;
		loadVector->mat[meshDb.elementDb[i].nodeId[1] - 1][0] += val;
		loadVector->mat[meshDb.elementDb[i].nodeId[2] - 1][0] += val;
	}
}
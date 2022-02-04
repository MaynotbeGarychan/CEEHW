/************************************************************
FileName: femMain.c
Author: Chen
Date: 2022/02/02
Description: main program of the fem analysis, 
			 the whole route was assemble in this src file
***********************************************************/
#pragma once
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "mesh.h"
#include "matrix.h"
#include "solver.h"
#include "Io.h"
#include "test.h"
#include "translator.h"
#include "model.h"
#include "femMain.h"

int femMain(Io ioInfo)
{
	// read txt
	mesh meshDb;
	analysis analysisInfo;
	readTxt(ioInfo, &meshDb, &analysisInfo);

	// translate information, set mesh info and solver param
	translator(&meshDb, &analysisInfo);

	// make a resultdb to get result
	result resultDb;

	// go to the specific problem
	switch (analysisInfo.solveProblem)
	{
	case WAVE_PRO:
		if (analysisInfo.usedTimeInteScheme)
		{

		}
		else
		{
			fem1dWaveStatic(ioInfo, meshDb, analysisInfo, &resultDb);
		}
		break;
	case POIS_PRO:
		fem2dPoissonStatic(ioInfo, meshDb, analysisInfo, &resultDb);
	default:
		break;
	}

	// write results to txt
	writeTxt(ioInfo, meshDb, analysisInfo, resultDb);

	return 1;
}

int fem1dWaveStatic(Io ioInfo, mesh meshDb, analysis analysisInfo, result* resultDb)
{

}

int fem2dPoissonStatic(Io ioInfo, mesh meshDb, analysis analysisInfo, result *resultDb)
{
	// for example, if we want to use newmark beta, we can set its step, 
	// but if it's simple fem problem, we can set the step is 1
	// assemble stiffness matrix
	matrix linearSystem;
	initilizeMatrix(&linearSystem, meshDb.meshInfoDb.nodeNum, meshDb.meshInfoDb.nodeNum + 1);
	assemble2DPoissonStatic(meshDb, analysisInfo, &linearSystem);

	// applied boundary conditions
	matrixInt idArray;
	allocateMatrixInt(&idArray, meshDb.meshInfoDb.nodeNum, 1);
	initializeIdArray(meshDb, &idArray);

	matrix slimLinearSys;
	allocateMatrix(&slimLinearSys, analysisInfo.dof[0], meshDb.meshInfoDb.nodeNum + 1);
	matrixInt slimIdArray;
	allocateMatrixInt(&slimIdArray, analysisInfo.dof[0], 1);
	deleteBoundaryRows(meshDb.meshInfoDb.nodeNum, analysisInfo.dof[0], analysisInfo.internalNodeIdList[0],
		linearSystem, idArray, &slimLinearSys, &slimIdArray);
	//freeMatrix(&linearSystem);
	//freeMatrixInt(&idArray);

	matrix finalLinearSys;
	allocateMatrix(&finalLinearSys, analysisInfo.dof[0], analysisInfo.dof[0] + 1);
	applyBoundaryConditionAndDeleteCols(meshDb.boundaryInfoDb.staticBoundaryNum, meshDb.boundaryInfoDb.dynamicBoundaryNumStep[0],
		analysisInfo.dof[0], meshDb.staticBoundaryDb, meshDb.dynamicBoundaryDb[0], analysisInfo.internalNodeIdList[0],
		slimLinearSys, &finalLinearSys);

	// call matrix solver
	matrix resultArr;
	allocateMatrix(&resultArr, analysisInfo.dof[0], 1);
	callMatrixSolver(analysisInfo.solverParam.matrixSolverType, finalLinearSys, slimIdArray, &resultArr);

	// save results
	saveScalarResultStep(0, 0, meshDb, analysisInfo, resultArr, slimIdArray, resultDb);

	return 1;
}
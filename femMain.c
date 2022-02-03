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
	allocateMatrix(&slimLinearSys, analysisInfo.dof, meshDb.meshInfoDb.nodeNum + 1);
	matrixInt slimIdArray;
	allocateMatrixInt(&slimIdArray, analysisInfo.dof, 1);
	deleteBoundaryRows(meshDb, analysisInfo, linearSystem, idArray, &slimLinearSys, &slimIdArray);
	//freeMatrix(&linearSystem);
	freeMatrix(&idArray);

	matrix finalLinearSys;
	allocateMatrix(&finalLinearSys, analysisInfo.dof, analysisInfo.dof + 1);
	applyBoundaryConditionAndDeleteCols(meshDb, analysisInfo, slimLinearSys, &finalLinearSys);

	// call matrix solver
	matrix resultArr;
	allocateMatrix(&resultArr, analysisInfo.dof, 1);
	switch (analysisInfo.solverParam.matrixSolverType)
	{
	case GAUSSPIVOT_SOLVER:
		gaussianEliminationSolveMatrix(finalLinearSys, &slimIdArray, &resultArr);
		break;
	case CG_SOLVER:
		conjugateGradientSolveMatrix(finalLinearSys, 1e-3, &resultArr);
		break;
	default:
		gaussianEliminationSolveMatrix(finalLinearSys, &slimIdArray, &resultArr);
		break;
	}

	// output results



	return 1;
}

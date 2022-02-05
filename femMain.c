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

	// make a resultDb to get result
	result resultDb;

	// go to the specific problem
	if (analysisInfo.usedTimeInteScheme == 1)
	{
		femMainDynamic(ioInfo, meshDb, analysisInfo, &resultDb);
	}
	else
	{
		femMainStatic(ioInfo, meshDb, analysisInfo, &resultDb);
	}

	// write results to txt
	writeTxt(ioInfo, meshDb, analysisInfo, resultDb);

	return 1;
}

int femMainStatic(Io ioInfo, mesh meshDb, analysis analysisInfo, result *resultDb)
{
	// assemble stiffness matrix
	matrix linearSystem;
	initilizeMatrix(&linearSystem, meshDb.meshInfoDb.nodeNum, meshDb.meshInfoDb.nodeNum + 1);

	switch (analysisInfo.solveProblem)
	{
	case POIS_PRO:
		assemble2DPoissonStatic(meshDb, analysisInfo, &linearSystem);
		break;
	case WAVE_PRO:
		assemble1DWaveStatic(meshDb, analysisInfo, &linearSystem);
		break;
	default:
		break;
	}

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

	matrix newDisplacementVec;
	allocateMatrix(&newDisplacementVec, meshDb.meshInfoDb.nodeNum, 1);
	assembleNewDisplacementVec(resultArr, slimIdArray, &newDisplacementVec,
		meshDb.boundaryInfoDb.staticBoundaryNum, 0, meshDb.staticBoundaryDb, meshDb.dynamicBoundaryDb[0]);

	// save results
	saveScalarResultStep(0, 0, meshDb.nodeDb,idArray, newDisplacementVec,resultDb);

	return 1;
}

int femMainDynamic(Io ioInfo, mesh meshDb, analysis analysisInfo, result* resultDb)
{
	matrix oldDisplacementVec, newDisplacementVec, oldVelocityVec, newVelocityVec, oldAccelerationVec, newAccelearationVec;
	matrix massMatrix, linearSystem;

	initilizeMatrix(&oldDisplacementVec, meshDb.meshInfoDb.nodeNum, 1);
	initilizeMatrix(&newDisplacementVec, meshDb.meshInfoDb.nodeNum, 1);
	initilizeMatrix(&oldVelocityVec, meshDb.meshInfoDb.nodeNum, 1);
	initilizeMatrix(&newVelocityVec, meshDb.meshInfoDb.nodeNum, 1);
	initilizeMatrix(&oldAccelerationVec, meshDb.meshInfoDb.nodeNum, 1);
	initilizeMatrix(&newAccelearationVec, meshDb.meshInfoDb.nodeNum, 1);
	initilizeMatrix(&massMatrix, meshDb.meshInfoDb.nodeNum, meshDb.meshInfoDb.nodeNum);
	initilizeMatrix(&linearSystem, meshDb.meshInfoDb.nodeNum, meshDb.meshInfoDb.nodeNum + 1);

	matrixInt idArray;
	allocateMatrixInt(&idArray, meshDb.meshInfoDb.nodeNum, 1);

	int currTime = 0;
	for (int step = 0; step < analysisInfo.timeInteParam.stepNum; step++)
	{
		currTime = step * analysisInfo.timeInteParam.stepLength + analysisInfo.timeInteParam.startTime;
		// assemble stiffness matrix
		zeroMatrix(&linearSystem);

		switch (analysisInfo.solveProblem)
		{
		case POIS_PRO:
			
			break;
		case WAVE_PRO:
			assembleGlobalStiffnessMatrix1DWaveDynamic(meshDb.meshInfoDb, meshDb.elementDb, meshDb.nodeDb, &linearSystem);
			break;
		default:
			break;
		}

		// add the velocity matrix
		addVelocityMatrixToGlobalMatrix1DWaveDynamic(10, oldVelocityVec, &linearSystem);

		// add the acceleration matrix
		matrix elemMassMat[MAX_NUM_ELEM];
		for (int i = 0; i < meshDb.meshInfoDb.elementNum; i++)
		{
			initilizeMatrix(&elemMassMat[i], meshDb.meshInfoDb.elemNodeNum, meshDb.meshInfoDb.elemNodeNum);
			assembleElementMassMatrix1DWaveDynamic(meshDb.elementDb[i], meshDb.nodeDb, &elemMassMat[i]);
		}
		assembleGlobalStiffnessMatrix(meshDb.meshInfoDb, meshDb.elementDb, elemMassMat, &massMatrix);
		addMassMatrixToGlobalMatrix1DWaveDynamic(oldAccelerationVec, massMatrix, &linearSystem);

		// applied boundary conditions
		initializeIdArray(meshDb, &idArray);

		matrix slimLinearSys;
		allocateMatrix(&slimLinearSys, analysisInfo.dof[step], meshDb.meshInfoDb.nodeNum + 1);
		matrixInt slimIdArray;
		allocateMatrixInt(&slimIdArray, analysisInfo.dof[step], 1);
		deleteBoundaryRows(meshDb.meshInfoDb.nodeNum, analysisInfo.dof[step], analysisInfo.internalNodeIdList[step],
			linearSystem, idArray, &slimLinearSys, &slimIdArray);

		matrix finalLinearSys;
		allocateMatrix(&finalLinearSys, analysisInfo.dof[step], analysisInfo.dof[step] + 1);
		applyBoundaryConditionAndDeleteCols(meshDb.boundaryInfoDb.staticBoundaryNum, meshDb.boundaryInfoDb.dynamicBoundaryNumStep[step],
			analysisInfo.dof[step], meshDb.staticBoundaryDb, meshDb.dynamicBoundaryDb[step], analysisInfo.internalNodeIdList[step],
			slimLinearSys, &finalLinearSys);

		// call matrix solver
		matrix resultArr;
		allocateMatrix(&resultArr, analysisInfo.dof[step], 1);
		callMatrixSolver(analysisInfo.solverParam.matrixSolverType, finalLinearSys, slimIdArray, &resultArr);

		// copy the result to the new displacement
		assembleNewDisplacementVec(resultArr, slimIdArray, &newDisplacementVec, 
			meshDb.boundaryInfoDb.staticBoundaryNum, meshDb.boundaryInfoDb.dynamicBoundaryNumStep[step], meshDb.staticBoundaryDb, meshDb.dynamicBoundaryDb[step]);
		
		// calculate new velocity vector, acceleration vector
		updateVelocityAndAccelerationVector(0.25, analysisInfo.timeInteParam.stepLength, 
			oldDisplacementVec, newDisplacementVec, oldVelocityVec, &newVelocityVec, oldAccelerationVec, &newAccelearationVec);

		// save results
		saveScalarResultStep(step,step*analysisInfo.timeInteParam.stepLength, meshDb.nodeDb, idArray, newDisplacementVec, resultDb);

		// copy the old vec to new vec
		copyMatrix(newDisplacementVec, &oldDisplacementVec);
		copyMatrix(newVelocityVec, &oldVelocityVec);
		copyMatrix(newAccelearationVec, &oldAccelerationVec);

	}
	return 1;
}
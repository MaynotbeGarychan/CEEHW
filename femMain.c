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
	assemble1DWaveStatic(meshDb, analysisInfo, &linearSystem);

	// call matrix solver


	// output results



	return 1;
}

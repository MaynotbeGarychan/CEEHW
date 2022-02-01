/*
*		Main program for Finite Element Analysis
*       Header file: test.h
*		Author: CHEN Jiawei, the University of Tokyo
*		Date:	2022/01/26
*/
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

int femMain(Io ioInfo)
{
	// read txt
	mesh meshDb;
	analysis analysisInfo;
	readTxt(ioInfo, &meshDb, &analysisInfo);

	// translate information, set mesh info and solver param


	// select the specific problem



	// for example, if we want to use newmark beta, we can set its step, 
	// but if it's simple fem problem, we can set the step is 1
	// assemble stiffness matrix


	// call matrix solver


	// output results



	return 1;
}
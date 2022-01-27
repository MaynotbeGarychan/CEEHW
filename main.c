/*
*		Main program of the solution, to link with each module
*		Author: CHEN Jiawei, the University of Tokyo
*		Date:	2022/01/26
*/
#pragma once
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "test.h"

int main(void)
{
	enum mode
	{
		FEM_MAIN_MODE,FEM_TEST_MODE,MATRIX_TEST_MODE,MESH_TEST_MODE
	};
    int testMode = FEM_TEST_MODE;

	switch (testMode)
	{
	case FEM_TEST_MODE:
		femTest();
		break;
	case MATRIX_TEST_MODE:
		matrixTest();
		break;
	case MESH_TEST_MODE:
		meshDelauneyTest();
		break;
	default:
		break;
	}

    return 1;
}
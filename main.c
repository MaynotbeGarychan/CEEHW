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
		FEM_TEST_MODE,MATRIX_TEST_MODE,MESH_TEST_MODE
	};
    int testMode = FEM_TEST_MODE;

	switch (testMode)
	{
	case FEM_TEST_MODE:
		femTest();
		break;
	case MATRIX_TEST_MODE:
		matrixTest();
	case MESH_TEST_MODE:
		meshDelauneyTest();
		break;
	}

    return 1;
}
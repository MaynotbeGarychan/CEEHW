/************************************************************
FileName: main.c
Author: Chen
Date: 2022/02/01
Description: main program of this solution, 
			 translate the input command line and 
			 link with each platform.
***********************************************************/
#pragma once
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "Io.h"
#include "test.h"

int main(int argc, char* argv[])
{
	/*	read the command line, get the directory, select the function
	*	argv[0]: Executable file
	*   argv[1]: Function (fem, tri)
	*	argv[2]: Directory of input file
	*	argv[3]: Directory of output file
	*/
	if (argc != 4 && argv[1] != "mat")
	{
		printf("Input of command line is wrong!");
		return 1;
	}
	for (int i = 0; i < argc; i++)
	{
		printf("%s ", argv[i]);
	}
	// get the directory of input and output
	Io ioInfo;
	ioInfo.inputDir = argv[2];
	ioInfo.outputDir = argv[3];
	// select the function
	int runMode = 0;
	if (strcmp(argv[1], "fem") == 0)
	{
		runMode = FEM_MAIN_MODE;
	}
	else if (strcmp(argv[1], "tri") == 0)
	{
		runMode = DELAUTRI_MAIN_MODE;
	}
	else if (strcmp(argv[1], "mat") == 0)
	{
		runMode = MATRIX_TEST_MODE;
	}
	else
	{
		printf("Current mode unsupported!");
		return 1;
	}

	/*	begin to link with each platform
	*	
	*/
	int retCode;
	switch (runMode)
	{
	case FEM_MAIN_MODE:
		retCode = femMain(ioInfo);
		break;
	case FEM_TEST_MODE:
		retCode = femTest(ioInfo);
		break;
	case MATRIX_TEST_MODE:
		retCode = matrixTest();
		break;
	case DELAUTRI_MAIN_MODE:
		retCode = meshDelauneyTest(ioInfo);
		break;
	default:
		break;
	}

    return 0;
}
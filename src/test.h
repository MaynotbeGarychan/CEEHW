#pragma once

enum mode
{
	FEM_MAIN_MODE, FEM_TEST_MODE, MATRIX_TEST_MODE, DELAUTRI_MAIN_MODE
};

//int femTest(Io ioInfo);
int matrixTest();
int meshDelauneyTest(Io ioInfo);
//int poissonMain(Io ioInfo);
//int femMain(Io ioInfo);
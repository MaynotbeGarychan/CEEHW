#pragma once
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "mesh.h"
#include "matrix.h"
#include "solver.h"


int femTest()
{
    //      reading mesh file     //
    // -------------------------  //
    struct meshInfo meshInfoDb;
    // please specify the input file here
    FILE* fileIo = fopen("report2Mesh2.txt", "rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    char readType[4];
    // reading overall information for mesh 
    fscanf(fileIo, "%s %d %d %d %d", readType, &(meshInfoDb.id), &(meshInfoDb.nodeNum), &(meshInfoDb.elementNum), &(meshInfoDb.boundaryNum));

    struct node nodeDb[1000];
    struct element elementDb[1000];
    struct boundary boundaryDb[1000];
    // read nodes information
    for (int i = 0; i < meshInfoDb.nodeNum; i++)
    {
        fscanf(fileIo, "%s %d %lf %lf", readType, &(nodeDb[i].id), &(nodeDb[i].x), &(nodeDb[i].y));
    }
    // read element information
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        fscanf(fileIo, "%s %d %d %d %d", readType, &(elementDb[i].id), &(elementDb[i].nodeId[0]), &(elementDb[i].nodeId[1]), &(elementDb[i].nodeId[2]));
    }
    // read boundary information
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        fscanf(fileIo, "%s %d %lf", readType, &(boundaryDb[i].nodeId), &(boundaryDb[i].value));
    }
    fclose(fileIo);

    // translate the information //
    // ------------------------- //
    int bounNodeIdVec[1000];
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        bounNodeIdVec[i] = boundaryDb[i].nodeId;
    }
    int dof = meshInfoDb.nodeNum - meshInfoDb.boundaryNum;
    int unknownNodeVec[1000];
    int unknownNodeVecLen = 0;
    for (int i = 0; i < meshInfoDb.nodeNum; i++)
    {
        if (search(nodeDb[i].id, bounNodeIdVec, meshInfoDb.boundaryNum) == 0)
        {
            unknownNodeVec[unknownNodeVecLen] = nodeDb[i].id;
            unknownNodeVecLen++;
        }
    }

    // assemble related matrix //
    // ----------------------- //
    // assemble element matrix
    matrix elemMatrix[1000];
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        allocateMatrix(&(elemMatrix[i]), 3, 3);
        assembleElementStiffnessMatrix(elementDb[i], nodeDb, &(elemMatrix[i]));
    }
    // assemble global stiffness matrix
    matrix globalMatrix;
    initilizeMatrix(&globalMatrix, meshInfoDb.nodeNum, meshInfoDb.nodeNum);
    assembleGlobalStiffnessMatrix(meshInfoDb, elementDb, elemMatrix, &globalMatrix);
    printMatrix(&globalMatrix);
    // assemble loadvector
    matrix loadVector;
    initilizeMatrix(&loadVector, meshInfoDb.nodeNum, 1);
#ifdef OptionProblem
    assembleLoadVector(meshInfoDb, elementDb, 6, &loadVector);
#endif
    // make a system
    matrix systemMatrix;
    newCombineMatrixCol(globalMatrix, loadVector, &systemMatrix);
    printMatrix(&systemMatrix);

    // begin to solve the matrix //
    // ------------------------- //
    // delect unused rows
    matrixInt idvec;
    allocateMatrixInt(&idvec, meshInfoDb.nodeNum, 1);
    for (int i = 0; i < idvec.numRow; i++)
    {
        idvec.mat[i][0] = nodeDb[i].id;
    }
    matrix unknownSystemMatrix;
    allocateMatrix(&unknownSystemMatrix, dof, systemMatrix.numCol);
    matrixInt unknownIdVec;
    allocateMatrixInt(&unknownIdVec, dof, 1);
    deleteBoundaryRows(meshInfoDb, boundaryDb, systemMatrix, idvec, &unknownSystemMatrix, &unknownIdVec);
    // apply boundary conditions
    matrix appliedSystemMatrix;
    applyBoundaryCondtion(unknownSystemMatrix, unknownNodeVec, meshInfoDb, boundaryDb, &appliedSystemMatrix);
    printMatrix(&appliedSystemMatrix);
    // solve
    gaussianEliminationFEM(&appliedSystemMatrix, &unknownIdVec);
    printMatrix(&appliedSystemMatrix);

    // Output the result //
    // ----------------- //
    FILE* OutputIo = fopen("output.txt", "w");
    if (OutputIo == NULL)
    {
        return 0;
    }
    for (int i = 0; i < appliedSystemMatrix.numRow; i++)
    {
        int nodeId = unknownIdVec.mat[i][0];
        fprintf(fileIo, "x%d %lf %lf %lf\n", nodeId, nodeDb[nodeId - 1].x, nodeDb[nodeId - 1].y, appliedSystemMatrix.mat[i][appliedSystemMatrix.numCol - 1]);
    }
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        int nodeId = boundaryDb[i].nodeId;
        fprintf(fileIo, "x%d %lf %lf %lf\n", nodeId, nodeDb[nodeId - 1].x, nodeDb[nodeId - 1].y, boundaryDb[i].value);
    }
    fclose(OutputIo);

    return 1;
}
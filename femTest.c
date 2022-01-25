#pragma once
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "mesh.h"
#include "matrix.h"
#include "solver.h"
#include "Io.h"

#define testIo


int femTest()
{
#ifndef testIo
    char* meshFile = "report2Mesh2.txt";
    struct meshInfo meshInfoDb;
    struct node nodeDb[100];
    struct element elementDb[100];
    struct boundary boundaryDb[100];
    readMesh(meshFile, &meshInfoDb, &nodeDb, &elementDb, &boundaryDb);
#endif

    // ----------------------------------------------------------  //
    //                      reading mesh file                      //
    // ----------------------------------------------------------  //
    FILE* fileIo = fopen("report2Mesh2.txt", "rt");
    if (fileIo == NULL)
    {
        return 0;
    }
    struct meshInfo meshInfoDb;
    readMeshInfo(fileIo, &meshInfoDb);
    struct node nodeDb[1000];
    for (int i = 0; i < meshInfoDb.nodeNum; i++)
    {
        readNode(fileIo, &nodeDb[i]);
    }
    struct element elementDb[1000];
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        readElem(fileIo, &elementDb[i]);
    }
    struct boundary boundaryDb[1000];
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        readBoundary(fileIo, &boundaryDb[i]);
    }
    fclose(fileIo);


    // ----------------------------------------------------------  //
    //                     Translate the info                      //
    // ----------------------------------------------------------  //
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

    // ----------------------------------------------------------  //
    //                  Assemble related matrix                    //
    // ----------------------------------------------------------  //
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
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        freeMatrix(&(elemMatrix[i]));
    }
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
    freeMatrix(&globalMatrix);
    freeMatrix(&loadVector);


    // ----------------------------------------------------------  //
    //                  Solve the matrix system                    //
    // ----------------------------------------------------------  //
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
    matrix result;
    allocateMatrix(&result, appliedSystemMatrix.numRow, 1);
    enum Solver
    {
        GAUSS, CG
    };
    int solverMethod = GAUSS;
    switch (solverMethod)
    {
    case GAUSS:
        gaussianEliminationFEM(&appliedSystemMatrix, &unknownIdVec, &result);
        printMatrix(&appliedSystemMatrix);
    case CG:
        conjugateSolveMatrix(appliedSystemMatrix, 1e-3, &result);
    default:
        break;
    }

    // ----------------------------------------------------------  //
    //                      Output the result                      //
    // ----------------------------------------------------------  //
    
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
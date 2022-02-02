/************************************************************
FileName: solver.c
Author: Chen
Date: 2022/02/02
Description: constitutive modelling, time integration
***********************************************************/
#pragma once
#include "solver.h"
#include "matrix.h"

void assembleElementStiffnessMatrix(struct element elementDb, struct node nodeDb[], matrix* elemMatrix)
{
    // calculate the transformation matrix
    matrix transMatrix;
    allocateMatrix(&transMatrix, 2, 2);
    assembleTransMatrix(elementDb, nodeDb, &transMatrix);
    //printMatrix(&transMatrix);
    // its related value
    double detJ = calculateDetMatrix22(transMatrix);
    inverseMatrix(&transMatrix);
    matrix transInvtransMat;
    transposeMatrix(transMatrix, &transInvtransMat);

    // element stiffness matrix component
    matrix dfaidr;
    allocateMatrix(&dfaidr, 3, 2);
    dfaidr.mat[0][0] = -1;
    dfaidr.mat[0][1] = -1;
    dfaidr.mat[1][0] = 1;
    dfaidr.mat[1][1] = 0;
    dfaidr.mat[2][0] = 0;
    dfaidr.mat[2][1] = 1;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            elemMatrix->mat[i][j] = 0.5 * detJ * funcElemMatrix(transInvtransMat, dfaidr, i, j);
        }
    }
    printMatrix(elemMatrix);
}

double funcElemMatrix(matrix tranInvJ, matrix dfaidr, int m, int n)
{
    double dr1dx = tranInvJ.mat[0][0];
    double dr2dx = tranInvJ.mat[0][1];
    double dr1dy = tranInvJ.mat[1][0];
    double dr2dy = tranInvJ.mat[1][1];

    double x1sady1s = tranInvJ.mat[0][0] * tranInvJ.mat[0][0] + tranInvJ.mat[0][1] * tranInvJ.mat[0][1];
    double x1x2ady1y2 = tranInvJ.mat[0][0] * tranInvJ.mat[0][1] + tranInvJ.mat[1][0] * tranInvJ.mat[1][1];
    double y1sady2s = tranInvJ.mat[1][0] * tranInvJ.mat[1][0] + tranInvJ.mat[1][1] * tranInvJ.mat[1][1];

    double retVal = x1sady1s * dfaidr.mat[m][0] * dfaidr.mat[n][0] +
        x1x2ady1y2 * (dfaidr.mat[m][0] * dfaidr.mat[n][1] + dfaidr.mat[m][1] * dfaidr.mat[n][0]) +
        y1sady2s * dfaidr.mat[m][1] * dfaidr.mat[n][1];

    return retVal;
}

void assembleTransMatrix(struct element elementDb, struct node nodeDb[], matrix* transMatrix)
{
    int nodeOnePos = elementDb.nodeId[0] - 1;
    int nodeTwoPos = elementDb.nodeId[1] - 1;
    int nodeThreePos = elementDb.nodeId[2] - 1;

    transMatrix->mat[0][0] = -nodeDb[nodeOnePos].x + nodeDb[nodeTwoPos].x;
    transMatrix->mat[0][1] = -nodeDb[nodeOnePos].x + nodeDb[nodeThreePos].x;
    transMatrix->mat[1][0] = -nodeDb[nodeOnePos].y + nodeDb[nodeTwoPos].y;
    transMatrix->mat[1][1] = -nodeDb[nodeOnePos].y + nodeDb[nodeThreePos].y;
}

int deleteBoundaryRows(struct meshInfo meshInfoDb, struct boundary boundaryDb[], matrix inputMat, matrixInt idVec,
    matrix* outMat, matrixInt* unknownIdVec)
{
    // boundary node id vector
    int bounNodeVec[1000];
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        bounNodeVec[i] = boundaryDb[i].nodeId;
    }
    int row = 0;
    for (int i = 0; i < inputMat.numRow; i++)
    {
        if (search(i + 1, bounNodeVec, meshInfoDb.boundaryNum) == 0)
        {
            for (int j = 0; j < outMat->numCol; j++)
            {
                outMat->mat[row][j] = inputMat.mat[i][j];
            }
            unknownIdVec->mat[row][0] = idVec.mat[i][0];
            row++;
            //printMatrix(&outMat);
        }
    }
    int dof = meshInfoDb.nodeNum - meshInfoDb.boundaryNum;
    if (dof != row)
    {
        return 0;
    }
    return 1;
}

void applyBoundaryCondtion(matrix inputMat, int unkownNodeIdVec[], struct meshInfo meshInfoDb, struct boundary boundaryDb[], matrix* outMat)
{
    int lvPos = inputMat.numCol - 1;
    for (int i = 0; i < meshInfoDb.boundaryNum; i++)
    {
        double val = boundaryDb[i].value;
        int pos = boundaryDb[i].nodeId - 1;
        for (int j = 0; j < inputMat.numRow; j++)
        {
            inputMat.mat[j][lvPos] -= inputMat.mat[j][pos] * val;
            inputMat.mat[j][pos] = 0.0;
        }
    }
    allocateMatrix(outMat, inputMat.numRow, inputMat.numRow + 1);
    for (int i = 0; i < outMat->numRow; i++)
    {
        outMat->mat[i][outMat->numCol - 1] = inputMat.mat[i][lvPos];
    }
    for (int i = 0; i < meshInfoDb.nodeNum - meshInfoDb.boundaryNum; i++)
    {
        for (int j = 0; j < outMat->numRow; j++)
        {
            int unknownNodePos = unkownNodeIdVec[i] - 1;
            outMat->mat[j][i] = inputMat.mat[j][unknownNodePos];
        }
    }
}

void assembleLoadVector(struct meshInfo meshInfoDb, struct element elementDb[], double RHSvalue, matrix* loadVector)
{
    double val = -RHSvalue / 6;
    initilizeMatrix(loadVector, meshInfoDb.nodeNum, 1);
    for (int i = 0; i < meshInfoDb.elementNum; i++)
    {
        loadVector->mat[elementDb[i].nodeId[0] - 1][0] += val;
        loadVector->mat[elementDb[i].nodeId[1] - 1][0] += val;
        loadVector->mat[elementDb[i].nodeId[2] - 1][0] += val;
    }
}

int search(int val, int vec[], int vecLen)
{
    for (int i = 0; i < vecLen; i++)
    {
        if (val == vec[i])
        {
            return 1;
        }
    }
    return 0;
}

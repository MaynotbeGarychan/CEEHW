#pragma once
#include <stdio.h>

struct node
{
    int id;
    double x;
    double y;
    //double z;
};

struct element
{
    int id;
    int nodeId[3];
};

struct boundary
{
    int nodeId;
    double value;
};

struct meshInfo
{
    int id;
    int nodeNum;
    int elementNum;
    int boundaryNum;
};
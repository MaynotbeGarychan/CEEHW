cmake_minimum_required(VERSION 2.4)

project(hw)

include_directories(${PROJECT_SOURCE_DIR})
add_executable(hw main.c DelaunayTriangulation.c femMain.c femMain.c Io.c main.c matrix.c matrixTest.c mesh.c model.c Poisson2DFEM.c solver.c translator.c)
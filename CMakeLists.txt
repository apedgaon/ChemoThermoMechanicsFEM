cmake_minimum_required(VERSION 3.10)

project(CoupledFieldFEM VERSION 1.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(HEADER_FILES CoupledFieldFEM.h FEMUtility.h)
add_executable(CoupledFieldFEM CoupledFieldProblem.cpp CoupledFieldFEM.cpp FEMUtility.cpp)
target_include_directories(CoupledFieldFEM PUBLIC
							"${EIGEN_INCLUDE}")
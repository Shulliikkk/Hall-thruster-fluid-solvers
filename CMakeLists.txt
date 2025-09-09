cmake_minimum_required(VERSION 3.10)  
project(fluid_models)        
file(COPY "magnetic_field_distribution.csv"
DESTINATION ${CMAKE_BINARY_DIR})

add_compile_options(-O3)
add_compile_options(-ffast-math)

add_executable(thruster
	src/params.cpp
	src/conduction.cpp
	src/discharge.cpp
	src/walls.cpp
	src/solve.cpp
	src/main.cpp) 

include_directories("include" "include/csv" "include/eigen" "include/JustInterp")

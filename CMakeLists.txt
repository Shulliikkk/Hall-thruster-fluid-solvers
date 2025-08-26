cmake_minimum_required(VERSION 3.10)  # Минимальная версия CMake
project(fluid_models)                 # Название проекта

file(COPY "magnetic_field_distribution.csv"
DESTINATION ${CMAKE_BINARY_DIR})

add_executable(thruster
	src/params.cpp
	src/conduction.cpp
	src/discharge.cpp
	src/walls.cpp
	src/solve.cpp
	src/main.cpp)  # Создаём исполняемый файл

include_directories("include" "include/csv" "include/eigen" "include/JustInterp") # "include/csv")  # Линкуем Eigen


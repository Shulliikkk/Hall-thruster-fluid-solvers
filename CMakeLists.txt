cmake_minimum_required(VERSION 3.10)  # Минимальная версия CMake
project(fluid_models)                 # Название проекта

file(COPY "ionization_rate_data.csv" 
	"neutral_rate_data.csv" 
DESTINATION ${CMAKE_BINARY_DIR})

add_executable(thruster
	src/electricity.cpp
	src/walls.cpp
	src/solver.cpp
	src/main.cpp)  # Создаём исполняемый файл

include_directories("include" "include/eigen" "include/interpolator" "include/csv")  # Линкуем Eigen


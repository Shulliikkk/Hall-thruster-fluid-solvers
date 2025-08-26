cmake_minimum_required(VERSION 3.25)  # Минимальная версия CMake
project(fluid_models)                 # Название проекта

find_package(Eigen3 REQUIRED)

file(COPY "magnetic_field_distribution.csv"
DESTINATION ${CMAKE_BINARY_DIR})

add_executable(thruster
	src/params.cpp
	src/conduction.cpp
	src/discharge.cpp
	src/walls.cpp
	src/solve.cpp
	src/main.cpp)  # Создаём исполняемый файл

target_include_directories(thruster PRIVATE "include" "include/csv" "include/eigen" "include/JustInterp")
# Линкуем Eigen
# https://libeigen.gitlab.io/eigen/docs-nightly/TopicCMakeGuide.html
target_link_libraries(thruster PRIVATE Eigen3::Eigen)


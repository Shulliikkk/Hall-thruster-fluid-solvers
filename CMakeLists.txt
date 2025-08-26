cmake_minimum_required(VERSION 3.25)  # Минимальная версия CMake
project(fluid_models)                 # Название проекта

# https://cmake.org/cmake/help/latest/module/FetchContent.html
include(FetchContent)
FetchContent_Declare(
	eigen
 	# ссылка для скачивания (можно git репозиторий, а не архив, архив быстрее)
	URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
	# аргументы для find_package, чтобы постараться сначала найти, а потом уже качать (нужен cmake 3.24 и выше)
	FIND_PACKAGE_ARGS NAMES Eigen3
)
# после этой команды Eigen скачается - проверь папку build/_deps
FetchContent_MakeAvailable(eigen)

file(COPY "magnetic_field_distribution.csv"
DESTINATION ${CMAKE_BINARY_DIR})

add_executable(thruster
	src/params.cpp
	src/conduction.cpp
	src/discharge.cpp
	src/walls.cpp
	src/solve.cpp
	src/main.cpp)  # Создаём исполняемый файл

target_include_directories(thruster PRIVATE "include" "include/csv" "include/JustInterp")
# Линкуем Eigen
# https://libeigen.gitlab.io/eigen/docs-nightly/TopicCMakeGuide.html
target_link_libraries(thruster PRIVATE Eigen3::Eigen)


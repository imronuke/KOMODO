cmake_minimum_required(VERSION 3.15)

project(komodo)
enable_language(Fortran)

file(GLOB SRC_sc ${PROJECT_SOURCE_DIR}/src/*.f90 )
add_executable(komodo ${SRC_sc})
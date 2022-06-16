cmake_minimum_required(VERSION 3.15)

# set(CMAKE_Fortran_COMPILER nvfortran )
set(CMAKE_Fortran_COMPILER gfortran )
project(komodo Fortran)

option(DEBUG "Debug Mode                                 "  off)

if(CMAKE_Fortran_COMPILER_ID  STREQUAL "GNU")
   set(f90flags "-cpp -fbacktrace -fPIC -ffree-line-length-none")

   if(DEBUG)
      set(CMAKE_Fortran_FLAGS "-g  -Wall -pedantic -fbounds-check -ffpe-trap=invalid,overflow,underflow ${f90flags}")
   else()
      set(CMAKE_Fortran_FLAGS "-O3 ${f90flags}")
   endif()

   add_definitions(-D__GNU)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
   set(f90flags "-fpp -std08 -assume byterecl -traceback -diag-disable 6009 -fPIC -list-line-length-none")

   if(DEBUG)
      set(CMAKE_Fortran_FLAGS "-g -warn -ftrapuv -fp-stack-check -check all -check:noarg_temp_created -fpe0 ${f90flags}")
   else()
      set(CMAKE_Fortran_FLAGS "-O3 ${f90flags}")
   endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
   set(f90flags "-Mfree -Mpreprocess -Mbackslash")

   if(DEBUG)
      set(CMAKE_Fortran_FLAGS "-traceback -Wall -Kieee -Ktrap=fp ${f90flags}")
   else()
      set(CMAKE_Fortran_FLAGS "-O3 -fast ${f90flags}")
   endif()

   add_definitions(-D__NVHPC)
endif()

file(GLOB SRC_sc ${PROJECT_SOURCE_DIR}/src/*.f90 )
add_executable(komodo ${SRC_sc})

# Check git and get commit id and branch name
if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
    find_package(Git)
    if(GIT_FOUND)
      # Get the current working branch
      execute_process(
         COMMAND git rev-parse --abbrev-ref HEAD
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         OUTPUT_VARIABLE GIT_BRANCH
         OUTPUT_STRIP_TRAILING_WHITESPACE)

		  # Get the head commit date (latest commit date)
		  execute_process(
		  	COMMAND git log -1 --format=%cd
		  	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
		  	OUTPUT_VARIABLE GIT_DATE
		  	OUTPUT_STRIP_TRAILING_WHITESPACE)
       
      # Get the latest commit hash
      execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    	
    	# Prints message
    	message(STATUS "Git commit hash : ${GIT_COMMIT_HASH}")
		  message(STATUS "Git commit time : ${GIT_DATE}")
		  message(STATUS "Git branch name : ${GIT_BRANCH}")

    	add_definitions(-D__GIT_COMMIT_HASH="${GIT_COMMIT_HASH}")
		  add_definitions(-D__GIT_DATE="${GIT_DATE}")
		  add_definitions(-D__GIT_BRANCH="${GIT_BRANCH}")

		#add preprocess definition
		add_definitions(-D__GIT)
    endif()
endif()
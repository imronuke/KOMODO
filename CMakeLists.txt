cmake_minimum_required(VERSION 3.15)

option(GPU    "Use OpenACC GPU          "  OFF)
option(DEBUG  "DEBUGGING MODE           "  OFF)

# Set for debug only
if(GPU)
   set(CMAKE_Fortran_FLAGS "-fast -ta=tesla:lineinfo -Minfo=accel,ccff -Mfree -Mpreprocess -Mbackslash")
   set(CMAKE_Fortran_COMPILER nvfortran )
else()
   if(DEBUG)
      set(CMAKE_Fortran_FLAGS "-cpp -fbacktrace -g -pedantic -fbounds-check -ffpe-trap=invalid,overflow,underflow")
      set(CMAKE_Fortran_COMPILER gfortran )
   else()
      set(CMAKE_Fortran_FLAGS "-cpp -O4")
      set(CMAKE_Fortran_COMPILER gfortran )
   endif()
endif()

project(komodo Fortran)

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
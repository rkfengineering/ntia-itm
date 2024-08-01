cmake_minimum_required(VERSION 3.25.1 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

project(ntia-itm)

# Require at least GCC v12 or better
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "12.0")
    message(FATAL_ERROR "Insufficient gcc version (< 12.0)")
  endif()
  message(STATUS "GNU Compiler Version = ${CMAKE_CXX_COMPILER_VERSION}")
endif()

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)

# Assuming GCC toolchain
add_compile_options(-Wall -Wextra -Wpedantic)

find_package(GTest REQUIRED)

# On default, assume that no one wants to build tests, example script, or install library locally
option(P452_BUILD_TESTS "Indicates whether unit tests for P452 should be built" OFF)
option(P452_COMPILE_COVERAGE "Indicates whether P452 should be compiled with code coverage" OFF)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  	if (NTIA_ITM_BUILD_TESTS)
		# Used in unit tests
		enable_testing()
		if (NTIA_ITM_COMPILE_COVERAGE)
			add_compile_options(--coverage)
			add_link_options(--coverage)  
    	endif()
  	endif()
endif() 

# Only build these targets if they haven't been already (important when including this project as a submodule)
if (NOT TARGET ITMLib)
  	add_subdirectory(ITM)
endif()

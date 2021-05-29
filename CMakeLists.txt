cmake_minimum_required(VERSION 3.10.2)

project(IBSLib VERSION 0.1.0)

set(CMAKE_BUILD_TYPE DEBUG)

set(GCC_OPENMP_COMPILE_FLAGS "-fopenmp -D use_openmp -Ofast -march=native -ffast-math")
set(GCC_OPENMP_LINK_FLAGS "-fopenmp")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_OPENMP_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_OPENMP_LINK_FLAGS}")

add_executable(${PROJECT_NAME} libtest.cpp)

# numerics lib
add_library(numfunc SHARED
    NumericFunctions/NumericFunctions
)

# coulomb log functions lib
add_library(clog SHARED
    CoulombLog/CoulombLogFunctions
)

target_include_directories(numfunc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/NumericFunctions)
target_include_directories(clog PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/CoulombLog)

target_link_libraries(clog numfunc)
target_link_libraries(${PROJECT_NAME} numfunc)
target_link_libraries(${PROJECT_NAME} clog)
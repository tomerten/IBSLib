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

# integrators lib
add_library(integrators SHARED
    Integrators/Integrators
)

# add models lib
add_library(models SHARED
    Models/Models
)
# add radiation damping lib
add_library(raddamp SHARED 
    RadiationDamping/RadiationDamping
)

# add twiss read lib
add_library(twiss SHARED Twiss/Twiss)

target_include_directories(numfunc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/NumericFunctions)
target_include_directories(clog PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/CoulombLog)
target_include_directories(integrators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Integrators)
target_include_directories(models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Models)
target_include_directories(twiss PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Twiss)
target_include_directories(raddamp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/RadiationDamping)

target_link_libraries(clog numfunc)
target_link_libraries(integrators clog numfunc)
target_link_libraries(raddamp numfunc)
target_link_libraries(models integrators clog numfunc)
#target_link_libraries(integrators clog)

target_link_libraries(${PROJECT_NAME} numfunc)
target_link_libraries(${PROJECT_NAME} clog)
target_link_libraries(${PROJECT_NAME} integrators)
target_link_libraries(${PROJECT_NAME} models)
target_link_libraries(${PROJECT_NAME} twiss)
target_link_libraries(${PROJECT_NAME} raddamp)
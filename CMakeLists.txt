cmake_minimum_required(VERSION 3.10.2)

project(IBSLib VERSION 0.1.0)

set(CMAKE_BUILD_TYPE DEBUG)

set(GCC_OPENMP_COMPILE_FLAGS "-fopenmp -D use_openmp -Ofast -march=native -ffast-math")
set(GCC_OPENMP_LINK_FLAGS "-fopenmp ")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_OPENMP_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_OPENMP_LINK_FLAGS}")

# adding pybind11
SET(MY_PYBIND ${CMAKE_CURRENT_SOURCE_DIR}/third_party/pybind11-2.5.0)
add_subdirectory(${MY_PYBIND})

add_executable(${PROJECT_NAME} libtest.cpp)
add_executable(DemoODE  libtest.cpp)
add_executable(DemoNumericFunctions  ${CMAKE_CURRENT_SOURCE_DIR}/Demos/DemoNumericFunctions.cpp)
add_executable(DemoRadDamping  ${CMAKE_CURRENT_SOURCE_DIR}/Demos/DemoRadDamping.cpp)
add_executable(DemoCoulombLog  ${CMAKE_CURRENT_SOURCE_DIR}/Demos/DemoCoulombLog.cpp)
add_executable(DemoIntegrators  ${CMAKE_CURRENT_SOURCE_DIR}/Demos/DemoIntegrators.cpp)
add_executable(DemoIBS  ${CMAKE_CURRENT_SOURCE_DIR}/Demos/DemoIBS.cpp)

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

# add ODE lib
add_library(ode SHARED OrdDiffEq/OrdDiffEq)

target_include_directories(numfunc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/NumericFunctions)
target_include_directories(clog PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/CoulombLog)
target_include_directories(integrators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Integrators)
target_include_directories(models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Models)
target_include_directories(twiss PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Twiss)
target_include_directories(raddamp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/RadiationDamping)
target_include_directories(ode PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/OrdDiffEq)

target_link_libraries(raddamp numfunc)
target_link_libraries(clog raddamp numfunc)
target_link_libraries(integrators clog numfunc)
target_link_libraries(models integrators clog numfunc)
target_link_libraries(ode models integrators clog raddamp numfunc twiss)

target_link_libraries(${PROJECT_NAME} numfunc)
target_link_libraries(${PROJECT_NAME} clog)
target_link_libraries(${PROJECT_NAME} integrators)
target_link_libraries(${PROJECT_NAME} models)
target_link_libraries(${PROJECT_NAME} twiss)
target_link_libraries(${PROJECT_NAME} raddamp)
target_link_libraries(${PROJECT_NAME} ode)

target_link_libraries(DemoNumericFunctions numfunc twiss raddamp)
target_link_libraries(DemoRadDamping numfunc twiss raddamp)
target_link_libraries(DemoCoulombLog clog numfunc twiss raddamp)
target_link_libraries(DemoIntegrators integrators clog numfunc twiss raddamp)
target_link_libraries(DemoIBS models integrators clog numfunc twiss raddamp)
target_link_libraries(DemoODE models ode integrators clog numfunc twiss raddamp)

pybind11_add_module(ibslib_pb 
    ibslib_pb.cpp
    SHARED Twiss/Twiss.cpp
    )
#target_link_libraries(ibslib_pb ${PROJECT_NAME})
cmake_minimum_required(VERSION 3.11)

project(molecular-dynamics-simulation)

# set(CMAKE_BUILD_TYPE "Release")
# set(EIGEN3_INCLUDE_DIR "external/Eigen3/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)
# add_compile_options(-ansi -pedantic -Wall -Wextra -Weffc++)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)
include_directories(external)
#target_link_libraries(my_mpi_target ${MPI_C_LIBRARIES})

set(COMMON_HDRS
        ${CMAKE_CURRENT_SOURCE_DIR}/header/neighbors.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/gupta.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/atoms.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/xyz.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/verlet.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/kineticAndTemp.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/ljDS.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/lj.h
        ${CMAKE_CURRENT_SOURCE_DIR}/header/bThermostat.h
)

set(COMMON_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/src/neighbors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gupta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/xyz.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/verlet.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kineticAndTemp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ljDS.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/lj.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/bThermostat.cpp
)

set(MYPROJECT
        ${CMAKE_CURRENT_SOURCE_DIR}/header/milestones_serial.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/milestones_serial.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/main_serial.cpp
)

add_subdirectory(tests)

if(MPI_BUILD) # set -DMPI_BUILD=True if running on UniCluster.
 
        find_package(MPI REQUIRED)
        include_directories(SYSTEM ${MPI_INCLUDE_PATH})

        set(MYPROJECT_HDRS_SRCS_PARALLEL
                ${CMAKE_CURRENT_SOURCE_DIR}/header/mpi_support.h
                ${CMAKE_CURRENT_SOURCE_DIR}/header/domain.h
                ${CMAKE_CURRENT_SOURCE_DIR}/header/milestones_parallel.h
                ${CMAKE_CURRENT_SOURCE_DIR}/src/domain.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/milestones_parallel.cpp

        )

        set(MYPROJECT_PARALLEL
                ${CMAKE_CURRENT_SOURCE_DIR}/main_parallel.cpp
        )
        add_executable(parallel_sim ${COMMON_SRCS} ${COMMON_HDRS} ${MYPROJECT_HDRS_SRCS_PARALLEL} ${MYPROJECT_PARALLEL})
        target_link_libraries(parallel_sim PUBLIC Eigen3::Eigen ${MPI_C_LIBRARIES} MPI::MPI_CXX)

else()

        add_executable(sim ${COMMON_SRCS} ${COMMON_HDRS} ${MYPROJECT})
        target_link_libraries(sim PUBLIC Eigen3::Eigen)

endif(MPI_BUILD)
unset(MPI_BUILD CACHE)

#add_compile_options(-lmpi_cxx)

#add_executable(myproject ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${MYPROJECT_executable_SRCS})
# add_executable(ih ${IH_HDRS} ${IH})
#target_link_libraries(myproject PUBLIC Eigen3::Eigen)

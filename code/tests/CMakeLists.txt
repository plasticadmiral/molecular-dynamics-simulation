include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(myproject_tests_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/test_ljDS.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test_bThermostat.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test_neighbors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test_gupta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test_verlet.cpp
)

include_directories(executable_simulation_tests ${CMAKE_SOURCE_DIR})
add_executable(executable_simulation_tests ${COMMON_SRCS} ${COMMON_HDRS} ${myproject_tests_SRCS})
target_link_libraries(executable_simulation_tests PUBLIC Eigen3::Eigen gtest gtest_main)

include(GoogleTest)
gtest_discover_tests(executable_simulation_tests)

add_test(NAME executable_simulation_tests COMMAND executable_simulation_tests)

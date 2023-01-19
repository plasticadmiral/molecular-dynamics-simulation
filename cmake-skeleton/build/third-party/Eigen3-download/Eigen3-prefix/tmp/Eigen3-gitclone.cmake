
if(NOT "/mnt/d/documents/molecular" IS_NEWER_THAN "dynamics/cmake-skeleton/build/external")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: 'dynamics/cmake-skeleton/build/external'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/mnt/d/documents/molecular"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/mnt/d/documents/molecular'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "dynamics/cmake-skeleton/build/external/Eigen3"  clone --no-checkout --origin "3.4.0" "/usr/bin/git" ""
    WORKING_DIRECTORY "Eigen3"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: '/usr/bin/git'")
endif()

execute_process(
  COMMAND "dynamics/cmake-skeleton/build/external/Eigen3"  checkout https://gitlab.com/libeigen/eigen.git --
  WORKING_DIRECTORY "Eigen3/"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'https://gitlab.com/libeigen/eigen.git'")
endif()

set(init_submodules origin)
if(init_submodules)
  execute_process(
    COMMAND "dynamics/cmake-skeleton/build/external/Eigen3"  submodule update --recursive --init TRUE
    WORKING_DIRECTORY "Eigen3/"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'Eigen3/'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/mnt/d/documents/molecular"
    "dynamics/cmake-skeleton/build/external"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: 'dynamics/cmake-skeleton/build/external'")
endif()


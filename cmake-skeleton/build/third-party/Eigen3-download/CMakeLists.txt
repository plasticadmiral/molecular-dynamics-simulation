
cmake_minimum_required(VERSION 3.1)
project(Eigen3-download NONE)
include(ExternalProject)
ExternalProject_Add(Eigen3
	SOURCE_DIR /mnt/d/documents/molecular dynamics/cmake-skeleton/build/external/Eigen3
	BINARY_DIR /mnt/d/documents/molecular dynamics/cmake-skeleton/build/third-party/Eigen3-download
	GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
	GIT_TAG 3.4.0
	CONFIGURE_COMMAND ""
	BUILD_COMMAND     ""
	INSTALL_COMMAND   ""
	TEST_COMMAND      ""
	)

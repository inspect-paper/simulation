cmake_minimum_required (VERSION 2.8)
project (sim-inspect)

# required for DCE
set(CMAKE_INSTALL_PREFIX "/usr")

# set(CMAKE_CXX_COMPILER "g++-4.9")


add_definitions (-Werror -Wall -g
		-O3
		-std=c++11
#		-DNDEBUG
#		-DMT19937
#		-DNC_DEBUG
		)
  


#SET (BOOST_ROOT "/usr/include/boost_1_59_0")
#SET (BOOST_INCLUDEDIR "/usr/include/boost_1_59_0/")
#SET (BOOST_LIBRARYDIR "/usr/include/boost_1_59_0/lib")
  
#find_package(Boost 1.59.0 COMPONENTS program_options REQUIRED)
#MESSAGE(STATUS "Variable is: " ${Boost_LIBRARY_DIRS})

#include_directories(${BOOST_INCLUDEDIR})
include_directories(time)
include_directories(.)

#link_directories(${Boost_LIBRARY_DIRS})

add_library(nc SHARED linear_algebra.cc network_coding.cc)

add_executable (tests unit_tests.cxx unit_tests_nc.cxx)

target_link_libraries(tests nc)


install(TARGETS nc LIBRARY DESTINATION lib)
install(FILES network_coding.hxx linear_algebra.hxx Matrix.h DESTINATION include/nc)

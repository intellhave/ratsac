cmake_minimum_required(VERSION 2.8)

project(usac)

find_package(LAPACK REQUIRED)
#find_package(Armadillo REQUIRED)


set(GLOG_LIBRARIES "/usr/lib/x86_64-linux-gnu/")
set(LP_SOLVE_LIB "/usr/lib/")

set(OSI_LINK_LIBS Osi Clp OsiClp CoinUtils)


##INCLUDE DIRECTORIES
include_directories(problems/)
include_directories(samples/)
include_directories(mcts/)
include_directories(toolbox/)
include_directories(quasiconvex/)

include_directories(./)
include_directories(/usr/include/glog/)
include_directories(/usr/inlcude/gflags/)
include_directories(/usr/include/)
include_directories(/usr/include/coin/)
#include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(/usr/include/eigen3/)

## LINK DIRECTORIES
link_directories(${LAPACK_LIBRARIES} ${GLOG_LIBRARIES} )

set(CMAKE_CXX_FLAGS "-std=c++14"  ${CMAKE_CXX_FLAGS})

file(GLOB SOURCES "*.cc"
  "problems/*.cc"
  "samples/common.cc"
  "mcts/*.cc"
  "toolbox/LinearFit.cc"
  "toolbox/HomographyFit.cc"
  "toolbox/Astar.cc"
  "toolbox/config.cc"
  "toolbox/extra.cc"
  "toolbox/quasiconvexFunctions.cc"
  "toolbox/QuasiconvexFit.cc"
  "toolbox/LinearProgramming.cc";  
  )



#add_executable(test_homog "samples/test_homog.cc" ${SOURCES})
#add_executable(test_fund "samples/test_fund.cc" ${SOURCES})
#add_executable(test_linear "samples/test_linear.cc" ${SOURCES})
add_executable(linear_experiment "experiments/linear_experiment.cc" ${SOURCES})
add_executable(homo_experiment "experiments/homography.cc" ${SOURCES})

#target_link_libraries(test_homog glog gflags config++  lapack ${ARMADILLO_LIBRARIES} ${OSI_LINK_LIBS} )
#target_link_libraries(test_fund glog gflags config++  lapack ${ARMADILLO_LIBRARIES})
#target_link_libraries(test_linear glog gflags config++ lapack ${ARMADILLO_LIBRARIES})
#target_link_libraries(linear_experiment glog gflags config++ lapack ${ARMADILLO_LIBRARIES} ${OSI_LINK_LIBS})
#target_link_libraries(homo_experiment glog gflags config++ lapack ${ARMADILLO_LIBRARIES} ${OSI_LINK_LIBS})

target_link_libraries(linear_experiment glog gflags config++ lapack ${OSI_LINK_LIBS})
target_link_libraries(homo_experiment glog gflags config++ lapack ${OSI_LINK_LIBS})

cmake_minimum_required(VERSION 3.22)

project(BiCG)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(main src/main.cpp src/BiCGSTAB.cpp)

include_directories(./include)

target_link_libraries(main Eigen3::Eigen)

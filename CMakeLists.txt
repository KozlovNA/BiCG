cmake_minimum_required(VERSION 3.22)

project(BiCG)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(main src/main.cpp)

include_directories(./include)

target_link_libraries(main Eigen3::Eigen)
# to change oprimization, run this command in shell
# cmake .. -DCMAKE_BUILD_TYPE=Release
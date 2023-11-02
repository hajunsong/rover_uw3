cmake_minimum_required(VERSION 3.0.2)

project(rover_uw3)

add_compile_options(-std=c++11)

find_package(Eigen3 REQUIRED)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)

file(GLOB SOURDCES src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURDCES})
# target_link_libraries(${PROJECT_NAME})
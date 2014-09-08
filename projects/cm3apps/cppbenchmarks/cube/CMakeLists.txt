include_directories(${CMAKE_SOURCE_DIR}/../dG3D/src ${PYTHON_INCLUDE_DIR})
add_executable(cube cube.cc)
target_link_libraries(cube dG3D)

include_directories(${CMAKE_SOURCE_DIR}/../dG3D/src ${PYTHON_INCLUDE_DIR})
add_definitions(-DNONLOCALGMSH)
add_executable(thermic cube.cc)
target_link_libraries(thermic dG3D)

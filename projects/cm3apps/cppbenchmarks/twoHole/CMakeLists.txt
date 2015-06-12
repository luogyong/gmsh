include_directories(${CMAKE_SOURCE_DIR}/../dG3D/src ${PYTHON_INCLUDE_DIR})
add_definitions(-DNONLOCALGMSH)
add_executable(twoHole twoHole.cc)
target_link_libraries(twoHole dG3D)

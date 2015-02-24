set(src
  Mapper.cpp
  MapEntity.cpp
  MapVertex.cpp
  MapElement.cpp
  GeoExtractor.cpp

  Mesh.cpp
  GroupOfElement.cpp

  GroupOfJacobian.cpp
  Jacobian.cpp
)

add_sources(geometry "${src}")

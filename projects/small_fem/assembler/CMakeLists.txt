set(src
  Dof.cpp

  SystemAbstract.cpp
  SystemEigen.cpp

  SystemPETSc.cpp
)

add_sources(assembler "${src}")

set(src
  FormulationPoisson.cpp
  FormulationNeumann.cpp
  FormulationSteadySlow.cpp
  FormulationPML.cpp

  FormulationLagrangeField.cpp
  FormulationFieldLagrange.cpp

  FormulationEMDA.cpp
  FormulationUpdateEMDA.cpp
  FormulationOO2.cpp
  FormulationUpdateOO2.cpp

  TermProjectionField.cpp
  TermProjectionGrad.cpp

  Quadrature.cpp
)

add_sources(formulation "${src}")

set(src
  FormulationMass.cpp
  FormulationStiffness.cpp
  FormulationProjection.cpp
  FormulationPoisson.cpp
  FormulationNeumann.cpp
  FormulationSteadySlow.cpp

  FormulationLagrangeField.cpp
  FormulationFieldLagrange.cpp

  FormulationEMDA.cpp
  FormulationUpdateEMDA.cpp
  FormulationOO2.cpp
  FormulationUpdateOO2.cpp

  TermFieldField.cpp
  TermGradGrad.cpp
  TermCurlCurl.cpp
  TermProjectionField.cpp
  TermProjectionGrad.cpp

  Quadrature.cpp
)

add_sources(formulation "${src}")

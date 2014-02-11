set(src
  FormulationPoisson.cpp
  FormulationProjectionVectorReal.cpp
  FormulationProjectionVectorComplex.cpp
  FormulationProjectionScalarReal.cpp
  FormulationProjectionScalarComplex.cpp

  FormulationLagrangeField.cpp
  FormulationFieldLagrange.cpp

  FormulationNeumann.cpp

  FormulationStiffness.cpp
  FormulationMass.cpp

  FormulationSteadySlow.cpp

  FormulationEMDA.cpp
  FormulationUpdateEMDA.cpp
  FormulationOO2.cpp
  FormulationUpdateOO2.cpp

  Term.cpp
  TermFieldField.cpp
  TermGradGrad.cpp
  TermCurlCurl.cpp
  TermProjectionField.cpp
  TermProjectionGrad.cpp

  Quadrature.cpp
)

add_sources(formulation "${src}")

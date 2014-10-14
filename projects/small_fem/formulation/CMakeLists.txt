set(src
  FormulationHelper.cpp

  FormulationPoisson.cpp
  FormulationSommerfeld.cpp
  FormulationSteadySlow.cpp
  FormulationPML.cpp

  FormulationLagrange.cpp
  FormulationLagrangeOne.cpp
  FormulationLagrangeTwo.cpp

  FormulationEMDA.cpp
  FormulationUpdateEMDA.cpp

  FormulationOO2.cpp
  FormulationUpdateOO2.cpp

  FormulationOSRC.cpp
  FormulationOSRCOne.cpp
  FormulationOSRCTwo.cpp
  FormulationOSRCThree.cpp
  FormulationOSRCFour.cpp
  FormulationUpdateOSRC.cpp

  FormulationOSRCScalar.cpp
  FormulationOSRCScalarOne.cpp
  FormulationOSRCScalarTwo.cpp
  FormulationOSRCScalarThree.cpp
  FormulationOSRCScalarFour.cpp
  FormulationUpdateOSRCScalar.cpp
  FormulationOSRCHelper.cpp

  FormulationJFLee.cpp
  FormulationJFLeeOne.cpp
  FormulationJFLeeTwo.cpp
  FormulationJFLeeThree.cpp
  FormulationJFLeeFour.cpp
  FormulationJFLeeFive.cpp
  FormulationJFLeeSix.cpp
  FormulationJFLeeSeven.cpp
  FormulationJFLeeEight.cpp
  FormulationUpdateJFLee.cpp

  TermProjectionField.cpp
  TermProjectionGrad.cpp

  Quadrature.cpp
)

add_sources(formulation "${src}")

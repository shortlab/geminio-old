/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ArcMaterial.h"
#include "MaterialParameters.h"

registerMooseObject("GeminioApp", ArcMaterial);

template<>
InputParameters validParams<ArcMaterial>()
{
  InputParameters params = validParams<Material>();

  // These parameters are simple constants from the input file
  params.addParam<Real>("ConstantTemperature", 298, "The temperature for a constant-T simulation, in K");

  // These parameters are used in diffusivity calculations
  params.addParam<Real>("ChromiumActivationEnergy", 0.33, "Activation energy for vacancy diffusion, in eV");
  params.addParam<Real>("ChromiumD0", 1e8, "Diffusivity prefactor for vacancies, in nm^2/s");

  // Many material properties scale with temperature (to be implemented later)
  // params.addRequiredCoupledVar("temperature", "Temperature of the CRUD for viscosity calculation (water temperature)");

  return params;
}

ArcMaterial::ArcMaterial(const InputParameters & parameters)
  : Material(parameters),
    _constant_temperature(getParam<Real>("ConstantTemperature")),
    _Cr_D0(getParam<Real>("ChromiumD0")),
    _Cr_activation_energy(getParam<Real>("ChromiumActivationEnergy")),
    _chromium_diffusivity_matprop(declareProperty<Real>("ChromiumDiffusivityMatProp")),
    _thermal_conductivity_matprop(declareProperty<Real>("ThermalConductivityMatProp"))
    // _T(coupledValue("temperature"))
{
}

void
ArcMaterial::computeQpProperties()
{
  // Comment to allow variable temperatures
  Real _T = _constant_temperature;  

  _chromium_diffusivity_matprop[_qp] = _Cr_D0 *
    std::exp(-_Cr_activation_energy / (MaterialParameters::kB * _T));

  _thermal_conductivity_matprop[_qp] = 10.0 - 0.001 * _T;
}

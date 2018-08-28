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

#include "ClusterDensity.h"

registerMooseObject("GeminioApp", ClusterDensity);

template<>
InputParameters validParams<ClusterDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_vars", "coupled vacancy type variables");
  params.addParam<Real>("scaling_factor", 1.0, "scaling factor to overall cluster density, eg. atomic volume");
  return params;
}

ClusterDensity::ClusterDensity(const InputParameters & parameters)
  : AuxKernel(parameters),
    _scaling_factor(getParam<Real>("scaling_factor"))
{
  const auto ncoupled = coupledComponents("coupled_vars");
  _no_vars.resize(ncoupled);
  _val_vars.resize(ncoupled);

  for (unsigned int i = 0; i < ncoupled; ++i)
  {
    _no_vars[i] = coupled("coupled_vars", i);
    _val_vars[i] = &coupledValue("coupled_vars", i);
  }
}

Real
ClusterDensity::computeValue()
{
  // total defect cluster concentration
  Real total = 0.0;
  const auto no_vars = _no_vars.size() ;
  for (unsigned int i = 0; i < no_vars; ++i)
    total += (*_val_vars[i])[_qp];

  return total * _scaling_factor;
}

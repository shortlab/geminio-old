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

#include "ConstantKernel.h"
#include <limits>

registerMooseObject("GeminioApp", ConstantKernel);

template<>
InputParameters validParams<ConstantKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("value", 0.0, "add kernel: -value * _test[_i][_qp]");
  params.addParam<Real>("tlimit", std::numeric_limits<Real>::max(), "set lifetime for the kernel");
  return params;
}

ConstantKernel::ConstantKernel(const InputParameters & parameters) :
    Kernel(parameters),
    _t_limit(getParam<Real>("tlimit")),
    _val(getParam<Real>("value"))
{
}

Real
ConstantKernel::computeQpResidual()
{
  if (_t < _t_limit)
    return -_val * _test[_i][_qp];
  return 0.0;
}

Real
ConstantKernel::computeQpJacobian()
{
  return 0.0;
}

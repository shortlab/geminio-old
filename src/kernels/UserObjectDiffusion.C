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

#include "UserObjectDiffusion.h"
#include "GeminioUtils.h"

registerMooseObject("GeminioApp", UserObjectDiffusion);

template<>
InputParameters validParams<UserObjectDiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addRequiredParam<UserObjectName>("user_object", "the name of user object providing group constant");
  params.addParam<Real>("coeff", 1.0, "coefficient");
  return params;
}

UserObjectDiffusion::UserObjectDiffusion(const
     InputParameters & parameters)
     :Diffusion(parameters),
     _coeff(getParam<Real>("coeff")),
     _gc(getUserObject<GroupConstant>("user_object"))
{
  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  groupNo = GeminioUtils::getGroupNumber(cur_var_name);
}

Real
UserObjectDiffusion::computeQpResidual()
{
  Real gc = _gc._diff(groupNo);
  return  _coeff * gc * Diffusion::computeQpResidual();
}

Real
UserObjectDiffusion::computeQpJacobian()
{
  Real gc = _gc._diff(groupNo);
  return  _coeff * gc * Diffusion::computeQpJacobian();
}

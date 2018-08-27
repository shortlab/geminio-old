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

#include "AddGSumSIAClusterDensity.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddGSumSIAClusterDensity, "add_aux_kernel");

template<>
InputParameters validParams<AddGSumSIAClusterDensity>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<AuxVariableName>("aux_var", "aux variable name to hold value");
  params.addRequiredParam<UserObjectName>("group_constant", "User object holding the grouping method parameterization");
  params.addParam<Real>("scale_factor", 1.0, "scale factor used in the kernel");
  params.addParam<unsigned int>("lower_bound","starting size to count, inclusive");
  params.addParam<unsigned int>("upper_bound","ending size to count, inclusive");
  return params;
}

AddGSumSIAClusterDensity::AddGSumSIAClusterDensity(const InputParameters & params) :
    Action(params)
{
}

// only emission of point defects of same type are considered
void
AddGSumSIAClusterDensity::act()
{
  const auto number_i = getParam<unsigned int>("number_i");
  const auto scale_factor = getParam<Real>("scale_factor");

  const auto aux_var = getParam<AuxVariableName>("aux_var");
  const auto group_constant_name = getParam<UserObjectName>("group_constant");

  std::vector<VariableName> coupled_i_vars;
  for (unsigned int cur_num = 1; cur_num <= number_i; ++cur_num)
  {
    coupled_i_vars.push_back(name() + "0i" + Moose::stringify(cur_num));
    coupled_i_vars.push_back(name() + "1i" + Moose::stringify(cur_num));
  }

  InputParameters params = _factory.getValidParams("GSumSIAClusterDensity");
  params.set<AuxVariableName>("variable") = aux_var;
  params.set<std::vector<VariableName> > ("coupled_vars") = coupled_i_vars;
  params.set<Real>("scale_factor") = scale_factor;
  if (isParamValid("lower_bound"))
    params.set<unsigned int>("lower_bound") = getParam<unsigned int>("lower_bound");
  if (isParamValid("upper_bound"))
    params.set<unsigned int>("upper_bound") = getParam<unsigned int>("upper_bound");
  params.set<UserObjectName>("user_object") = group_constant_name;
  _problem->addAuxKernel("GSumSIAClusterDensity", "GSumSIAClusterDensity_" + aux_var, params);
}

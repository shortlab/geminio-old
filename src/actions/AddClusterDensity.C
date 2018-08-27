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

#include "AddClusterDensity.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "Conversion.h"
#include "ClusterDensity.h"

registerMooseAction("GeminioApp", AddClusterDensity, "add_aux_kernel");

template<>
InputParameters validParams<AddClusterDensity>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<AuxVariableName>("aux_var", "aux variable name to hold value");
  params.addRequiredParam<std::string>("var_prefix", "The prefix string of variables (in front of number)");
  params.addRequiredParam<std::vector<unsigned int>>("size_range", "The number of vacancy variables to solve");
  params.addParam<Real>("scale_factor", 1.0, "A scale factor to be applied to the variable");
  return params;
}

AddClusterDensity::AddClusterDensity(const InputParameters & params) :
    Action(params),
    _size_range(getParam<std::vector<unsigned int> >("size_range"))
{
  if (_size_range.size() != 2)
    paramError("size_range", "A vector with two entries denoting the size range is required");
}

void
AddClusterDensity::act()
{
  const auto aux_var = getParam<AuxVariableName>("aux_var");
  const auto var_prefix = getParam<std::string>("var_prefix");
  const auto scale_factor = getParam<Real>("scale_factor");

  std::vector<VariableName> coupled_vars;
  for (unsigned int i = _size_range[0]; i <= _size_range[1]; ++i)
    coupled_vars.push_back(var_prefix + Moose::stringify(i));

  InputParameters params = _factory.getValidParams("ClusterDensity");
  params.set<AuxVariableName>("variable") = aux_var;
  params.set<std::vector<VariableName> > ("coupled_vars") = coupled_vars;
  params.set<Real>("scaling_factor") = scale_factor;
  _problem->addAuxKernel("ClusterDensity", "ClusterDensity_" + aux_var, params);
}

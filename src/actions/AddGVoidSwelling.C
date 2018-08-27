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

#include "AddGVoidSwelling.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddGVoidSwelling, "add_kernel");

template<>
InputParameters validParams<AddGVoidSwelling>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<AuxVariableName>("aux_var", "aux variable name to hold value");
  params.addRequiredParam<UserObjectName>("group_constant", "User object holding the grouping method parameterization");
  return params;
}

AddGVoidSwelling::AddGVoidSwelling(const InputParameters & params) :
    Action(params)
{
}

// only emission of point defects of same type are considered
void
AddGVoidSwelling::act()
{
  const auto number_v = getParam<unsigned int>("number_v");

  const auto aux_var = getParam<AuxVariableName>("aux_var");
  const auto group_constant_name = getParam<UserObjectName>("group_constant");

  std::vector<VariableName> coupled_v_vars;
  for (int cur_num = 1; cur_num <= number_v; cur_num++)
  {
    coupled_v_vars.push_back(name() +"0v" + Moose::stringify(cur_num));
    coupled_v_vars.push_back(name() +"1v" + Moose::stringify(cur_num));
  }

  InputParameters params = _factory.getValidParams("GVoidSwelling");
  params.set<AuxVariableName>("variable") = aux_var;
  params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
  params.set<UserObjectName>("user_object") = group_constant_name;
  _problem->addAuxKernel("GVoidSwelling", "GVoidSwelling_" + aux_var, params);
}

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

#include "AddGVariable.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddGVariable, "add_variable");
registerMooseAction("GeminioApp", AddGVariable, "add_bc");
registerMooseAction("GeminioApp", AddGVariable, "add_ic");

template<>
InputParameters validParams<AddGVariable>()
{
  InputParameters params = validParams<GeminioAddVariableAction>();
  params.addParam<std::vector<unsigned int> >("IC_v_size","vacancy species number with initial concentration not ZERO");
  params.addParam<std::vector<unsigned int> >("IC_i_size","interstitial species number with initial concentration not ZERO");
  params.addParam<std::vector<Real> >("IC_v", "initial value for vacancy cluster corresponding to IC_v_size");
  params.addParam<std::vector<Real> >("IC_i", "initial value for interstitial cluster corresponding to IC_i_size");
  return params;
}

AddGVariable::AddGVariable(const InputParameters & params) :
    GeminioAddVariableAction(params)
{
}

void
AddGVariable::act()
{
  const auto vv = getParam<std::vector<unsigned int> >("IC_v_size");
  const auto ii = getParam<std::vector<unsigned int> >("IC_i_size");
  const auto initial_v = getParam<std::vector<Real> >("IC_v");
  const auto initial_i = getParam<std::vector<Real> >("IC_i");
  if (vv.size() != initial_v.size() || ii.size() != initial_i.size())
    mooseError("IC_v_size and IC_v should have same length, so are IC_i_size and IC_i., groupsize = 1 ");

  if (_current_task == "add_variable")
  {
    for (int cur_num = 1; cur_num <= _number_v; cur_num++)
    {
      addVariable(name() + "0v" + Moose::stringify(cur_num));
      addVariable(name() + "1v" + Moose::stringify(cur_num));
    }

    for (int cur_num = 1; cur_num <= _number_i; cur_num++)
    {
      addVariable(name() + "0i" + Moose::stringify(cur_num));
      addVariable(name() + "1i" + Moose::stringify(cur_num));
    }
  }

  else if (_current_task == "add_bc")
  {
    std::string var_name;
    for (unsigned int cur_num = 1; cur_num <= _number_v; ++cur_num)
    {
      var_name = name() + "0v" + Moose::stringify(cur_num);
      InputParameters params = _factory.getValidParams(_bc_name);
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params.set<Real>("value") = _boundary_value;
      _problem->addBoundaryCondition(_bc_name, var_name + "_left", params);
      params.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params.set<Real>("value") = _boundary_value;
      _problem->addBoundaryCondition(_bc_name, var_name + "_right", params);

      var_name = name() + "1v" + Moose::stringify(cur_num);
      InputParameters params1 = _factory.getValidParams(_bc_name);
      params1.set<NonlinearVariableName>("variable") = var_name;
      params1.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params1.set<Real>("value") = _boundary_value;
      _problem->addBoundaryCondition(_bc_name, var_name + "_left", params1);
      params1.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params1.set<Real>("value") = _boundary_value;
      _problem->addBoundaryCondition(_bc_name, var_name + "_right", params1);
    }

    for (int cur_num = 1; cur_num <= _number_i; cur_num++)
    {
      var_name = name() + "0i" + Moose::stringify(cur_num);
      InputParameters params = _factory.getValidParams(_bc_name);
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params.set<Real>("value") = _boundary_value;
      _problem->addBoundaryCondition(_bc_name, var_name + "_left", params);
      params.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params.set<Real>("value") = _boundary_value;
      _problem->addBoundaryCondition(_bc_name, var_name + "_right", params);

      var_name = name() + "1i" + Moose::stringify(cur_num);
      InputParameters params1 = _factory.getValidParams(_bc_name);
      params1.set<NonlinearVariableName>("variable") = var_name;
      params1.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params1.set<Real>("value") = _boundary_value;
      _problem->addBoundaryCondition(_bc_name, var_name + "_left", params1);
      params1.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params1.set<Real>("value") = _boundary_value;
      _problem->addBoundaryCondition(_bc_name, var_name + "_right", params1);
    }
  }

  else if (_current_task == "add_ic")
  {
    std::string var_name;
    for (unsigned int cur_num = 1; cur_num <= _number_v; ++cur_num)
    {
      auto it = std::find(vv.begin(), vv.end(), cur_num);
      addConstantIC(name() + "0v" + Moose::stringify(cur_num), it == vv.end() ? 0.0 : initial_v[it - vv.begin()]);
      addConstantIC(name() + "1v" + Moose::stringify(cur_num), 0.0);
    }

    for (int cur_num = 1; cur_num <= _number_i; ++cur_num)
    {
      auto it = std::find(ii.begin(), ii.end(), cur_num);
      addConstantIC(name() + "0i" + Moose::stringify(cur_num), it == ii.end() ? 0.0 : initial_v[it - ii.begin()]);
      addConstantIC(name() + "1i" + Moose::stringify(cur_num), 0.0);
    }
  }

}

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

#include "AddLotsOfVariableAction.h"
#include "FEProblem.h"
#include "MooseUtils.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddLotsOfVariableAction, "add_variable");
registerMooseAction("GeminioApp", AddLotsOfVariableAction, "add_bc");
registerMooseAction("GeminioApp", AddLotsOfVariableAction, "add_ic");

template<>
InputParameters validParams<AddLotsOfVariableAction>()
{
  InputParameters params = validParams<GeminioAddVariableAction>();
  params.addParam<Real>("initial_condition", 0.0, "Specifies the initial condition for this variable");
  params.addParam<Real>("scaling", 1.0, "Specifies a scaling factor to apply to this variable");
  params.addParam<bool>("use_constIC", false, "Specifies whether use single value for all variables as initial condition");
  return params;
}

AddLotsOfVariableAction::AddLotsOfVariableAction(const InputParameters & params) :
    GeminioAddVariableAction(params)
{
}

void
AddLotsOfVariableAction::act()
{
  const auto use_constIC = getParam<bool>("use_constIC");

  if (_current_task == "add_variable")
  {
    // vacancy variables
    for (int cur_num = 1; cur_num <= _number_v; ++cur_num)
      addVariable(name() + "v" + Moose::stringify(cur_num));
    
    // interstitial variables
    for (int cur_num = 1; cur_num <= _number_i; ++cur_num)
      addVariable(name() + "i" + Moose::stringify(cur_num));
  }

  else if (_current_task == "add_bc")
  {
    for (unsigned int cur_num = 1; cur_num <= _number_v + _number_i; ++cur_num)
    {
      std::string var_name;
      if (cur_num <= _number_v)
        var_name = name() + "v" + Moose::stringify(cur_num);
      else
        var_name = name() + "i" + Moose::stringify(cur_num - _number_v);

      InputParameters params = _factory.getValidParams(_bc_name);
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<BoundaryName> >("boundary").push_back("left");
      params.set<Real>("value") = _boundary_value;

      _problem->addBoundaryCondition(_bc_name, var_name + "_left", params);

      params.set<std::vector<BoundaryName> >("boundary")[0] = "right";
      params.set<Real>("value") = _boundary_value;

      _problem->addBoundaryCondition(_bc_name, var_name + "_right", params);
    }
  }

  else if (use_constIC == true && _current_task == "add_ic")
  {
    Real initial = getParam<Real>("initial_condition");
    if (!MooseUtils::absoluteFuzzyEqual(initial, 0.0))
    {
      for (int cur_num = 1; cur_num <= _number_v; ++cur_num)
        addConstantIC(name() + "v" + Moose::stringify(cur_num), initial);
        
      for (int cur_num = 1; cur_num <= _number_i; ++cur_num)
        addConstantIC(name() + "i" + Moose::stringify(cur_num), initial);
    }
  }

}

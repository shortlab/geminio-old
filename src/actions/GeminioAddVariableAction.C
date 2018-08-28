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

#include "GeminioAddVariableAction.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "Conversion.h"

template<>
InputParameters validParams<GeminioAddVariableAction>()
{
  InputParameters params = validParams<AddVariableAction>();

  params.addRequiredParam<unsigned int>("number_v", "Total number of vacancy group to add");
  params.addRequiredParam<unsigned int>("number_i", "Total number of interstitial group to add");

  MooseEnum bc_type_enum("neumann=0 dirichlet", "neumann");
  params.addParam<MooseEnum>("bc_type", bc_type_enum, "Boundary conditions to set for the generated variables, depending on spatial dependence");

  params.addParam<Real>("boundary_value", 0.0, "Specifies the initial condition for this variable");
  return params;
}

GeminioAddVariableAction::GeminioAddVariableAction(const InputParameters & params) :
    AddVariableAction(params),
    _bc_type(getParam<MooseEnum>("bc_type")),
    _boundary_value(getParam<Real>("boundary_value")),
    _number_i(getParam<unsigned int>("number_i")),
    _number_v(getParam<unsigned int>("number_v"))
{
  switch (_bc_type)
  {
    case 0:
      _bc_name = "NeumannBC";
      break;

    case 1:
      _bc_name = "DirichletBC";
      break;

    default:
      paramError("bc_type", "Invalid boundary condition.");
  }

  if (_scalar_var)
    paramError("family", "This action does not support adding scalar variables");
}

void
GeminioAddVariableAction::addConstantIC(const std::string & var_name, Real initial)
{
  InputParameters params = _factory.getValidParams("ConstantIC");
  params.set<VariableName>("variable") = var_name;
  params.set<Real>("value") = initial;
  _problem->addInitialCondition("ConstantIC", var_name + "_ic", params);
}

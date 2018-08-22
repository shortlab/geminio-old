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

#include "AddUserObjectDislocationSink.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "UserObjectSingleVariable.h"

#include <sstream>
#include <stdexcept>
#include <algorithm>
// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

template<>
InputParameters validParams<AddUserObjectDislocationSink>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<unsigned int> >("mobile_v_size", "A vector of mobile species");
  params.addRequiredParam<std::vector<unsigned int> >("mobile_i_size", "A vector of mobile species");
  params.addRequiredParam<std::string>("group_constant", "user object name");
  params.addRequiredParam<std::string>("call_function", "map to function to call in the userobject");
  return params;
}

AddUserObjectDislocationSink::AddUserObjectDislocationSink(const InputParameters & params) :
    AddVariableAction(params)
{
}

//only emission of point defects of same type are considered
void
AddUserObjectDislocationSink::act()
{
  const auto number_v = getParam<unsigned int>("number_v");
  const auto number_i = getParam<unsigned int>("number_i");
  std::vector<unsigned int> v_size = getParam<std::vector<unsigned int> >("mobile_v_size");
  std::vector<unsigned int> i_size = getParam<std::vector<unsigned int> >("mobile_i_size");
  std::string uo = getParam<std::string>("group_constant");
  std::string call_function = getParam<std::string>("call_function");
  std::vector<Real> vv,ii;

  const auto total_mobile_v = v_size.size();
  for (unsigned int i = 0; i < total_mobile_v; ++i)
    //account for sign; sink "+"
    vv.push_back(1.0);

  const auto total_mobile_i = i_size.size();
  for (unsigned int i = 0; i < total_mobile_i; ++i)
    //account for sign; sink "+"
    ii.push_back(1.0);


  for (unsigned int cur_num = 0; cur_num < total_mobile_v; ++cur_num)
  {
    Real coef = vv[cur_num];//cur_num to vacancy,"+"
    // add kernel for mobile v
    std::string var_name_v = name() + "v" + Moose::stringify(v_size[cur_num]);
    InputParameters params = _factory.getValidParams("UserObjectSingleVariable");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<Real>("coeff") = coef; // loss "+"
    params.set<UserObjectName>("user_object") = uo;
    params.set<std::string>("call_function") = call_function;
    _problem->addKernel("UserObjectSingleVariable", "DislocationSink_" + var_name_v + Moose::stringify(cur_num), params);
    _console << "add UserObjectSingleVariable disl: " << var_name_v << ", coef: " << coef << '\n';
  }

  for (unsigned int cur_num = 0; cur_num < total_mobile_i; ++cur_num)
  {
    Real coef = ii[cur_num];
    // add kernel for mobile i
    std::string var_name_i = name() + "i" + Moose::stringify(i_size[cur_num]);
    InputParameters params = _factory.getValidParams("UserObjectSingleVariable");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<Real>("coeff") = coef; // loss "+"
    params.set<UserObjectName>("user_object") = uo;
    params.set<std::string>("call_function") = call_function;
    _problem->addKernel("UserObjectSingleVariable", "DislocationSink_" + var_name_i + Moose::stringify(cur_num), params);
    _console << "add UserObjectSingleVariable disl: " << var_name_i << ", coef: " << coef << '\n';
  }
}

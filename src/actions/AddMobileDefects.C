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

#include "AddMobileDefects.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddMobileDefects, "add_kernel");

template<>
InputParameters validParams<AddMobileDefects>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<unsigned int>("max_mobile_v", "maximum size of mobile vacancy cluster");
  params.addRequiredParam<unsigned int>("max_mobile_i", "maximum size of mobile intersitial cluster");
  params.addRequiredParam<UserObjectName>("group_constant", "User object holding the grouping method parameterization");
  params.addParam<Real>("dislocation",0.0,"dislocation density");
  params.addParam<std::vector<Real> >("disl_mobile_v", "A vector of dislocation bias for mobile species");
  params.addParam<std::vector<Real> >("disl_mobile_i", "A vector of dislocation bias for mobile species");
  return params;
}

AddMobileDefects::AddMobileDefects(const InputParameters & params) :
    Action(params)
{
}

// only emission of point defects of same type are considered
void
AddMobileDefects::act()
{
  const auto number_v = getParam<unsigned int>("number_v");
  const auto number_i = getParam<unsigned int>("number_i");
  const auto max_mobile_v = getParam<unsigned int>("max_mobile_v");
  const auto max_mobile_i = getParam<unsigned int>("max_mobile_i");

  std::vector<unsigned int> v_size,i_size;
  for (unsigned int i = 0; i < max_mobile_v; ++i)
    v_size.push_back(i+1);
  for (unsigned int i = 0; i < max_mobile_i; ++i)
    i_size.push_back(i+1);

  const auto disl = getParam<Real>("dislocation");
  const auto group_constant_name = getParam<UserObjectName>("group_constant");
  std::vector<Real> disl_v = getParam<std::vector<Real> >("disl_mobile_v");
  std::vector<Real> disl_i = getParam<std::vector<Real> >("disl_mobile_i");


  const auto num_mobile_v = v_size.size();
  const auto num_mobile_i = i_size.size();
  // const auto disl_v_size = disl_v.size();
  // const auto disl_i_size = disl_i.size();

  std::vector<VariableName> coupled_v_vars;
  std::vector<VariableName> coupled_i_vars;
  std::string _prefix = name();
  std::string var_name;
  for (unsigned int i = 1; i <= number_v; ++i)
  {
    var_name = _prefix + "v" + Moose::stringify(i);
    coupled_v_vars.push_back(var_name);
  }
  for (unsigned int i = 1; i <= number_i; ++i)
  {
    var_name = _prefix + "i" + Moose::stringify(i);
    coupled_i_vars.push_back(var_name);
  }

  // first add mobile v
  for (unsigned int cur_num = 1; cur_num <= num_mobile_v; ++cur_num)
  {
    // dislocation absorption factor
    // (disl_v_size>0)?disl_v[cur_num-1]:1.0;
    Real coef = 1.0;
    std::string var_name_v = name() + "v" + Moose::stringify(v_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("MobileDefects");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<Real>("dislocation") = disl;
    params.set<Real>("dislocation_factor") = coef;
    params.set<UserObjectName>("user_object") = group_constant_name;
    params.set<unsigned int>("number_v") = number_v;
    params.set<unsigned int>("number_i") = number_i;
    params.set<std::vector<unsigned int> >("mobile_v_size") = v_size;
    params.set<std::vector<unsigned int> >("mobile_i_size") = i_size;
    _problem->addKernel("MobileDefects", "MobileDefects_" + var_name_v + "_" + Moose::stringify(cur_num), params);
    // printf("add MobileDefects: %s \n",var_name_v.c_str());
  }

  // second add mobile i
  for (unsigned int cur_num = 1; cur_num <= num_mobile_i; ++cur_num)
  {
    // dislocation absorption factor
    // (disl_i_size>0)?disl_i[cur_num-1]:1.0;
    Real coef = 1.0;
    std::string var_name_i = name() + "i" + Moose::stringify(i_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("MobileDefects");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<Real>("dislocation") = disl;
    params.set<Real>("dislocation_factor") = coef;
    params.set<UserObjectName>("user_object") = group_constant_name;
    params.set<unsigned int>("number_v") = number_v;
    params.set<unsigned int>("number_i") = number_i;
    params.set<std::vector<unsigned int> >("mobile_v_size") = v_size;
    params.set<std::vector<unsigned int> >("mobile_i_size") = i_size;
    _problem->addKernel("MobileDefects", "MobileDefects_" + var_name_i + "_" + Moose::stringify(cur_num), params);
    // printf("add MobileDefects: %s \n",var_name_i.c_str());
  }
}

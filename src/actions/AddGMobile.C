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

#include "AddGMobile.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

template<>
InputParameters validParams<AddGMobile>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<unsigned int>("max_mobile_v", "maximum size of mobile vacancy cluster");
  params.addRequiredParam<unsigned int>("max_mobile_i", "maximum size of mobile interstitial cluster");
  params.addRequiredParam<UserObjectName>("group_constant", "User object holding the grouping method parameterization");
  return params;
}

AddGMobile::AddGMobile(const InputParameters & params) :
    Action(params)
{
}

//only emission of point defects of same type are considered
void
AddGMobile::act()
{
  const auto number_v = getParam<unsigned int>("number_v");
  const auto number_i = getParam<unsigned int>("number_i");
  const auto num_mobile_v = getParam<unsigned int>("max_mobile_v");
  const auto num_mobile_i = getParam<unsigned int>("max_mobile_i");

  const auto group_constant_name = getParam<UserObjectName>("group_constant");

  std::vector<VariableName> coupled_v_vars;
  for (int cur_num = 1; cur_num <= number_v; ++cur_num)
  {
    coupled_v_vars.push_back(name() + "0v" + Moose::stringify(cur_num));
    coupled_v_vars.push_back(name() + "1v" + Moose::stringify(cur_num));
  }

  std::vector<VariableName> coupled_i_vars;
  for (int cur_num = 1; cur_num <= number_i; ++cur_num)
  {
    coupled_i_vars.push_back(name() + "0i" + Moose::stringify(cur_num));
    coupled_i_vars.push_back(name() + "1i" + Moose::stringify(cur_num));
  }

  // first add mobile v
  for (int cur_num = 1; cur_num <= num_mobile_v; ++cur_num)
  {
    std::string var_name_v = name() + "0v" + Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("GMobile");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<UserObjectName>("user_object") = group_constant_name;
    params.set<unsigned int>("number_v") = number_v;
    params.set<unsigned int>("number_i") = number_i;
    params.set<unsigned int>("max_mobile_v") = num_mobile_v;
    params.set<unsigned int>("max_mobile_i") = num_mobile_i;
    _problem->addKernel("GMobile", "GMobile_" + var_name_v + "_" + Moose::stringify(cur_num), params);
    //printf("add GMobile: %s \n",var_name_v.c_str());

    // add pesudo kernel for L1 coefficient
    var_name_v = name() + "1v" + Moose::stringify(cur_num);
    InputParameters params1 = _factory.getValidParams("ConstantKernel");
    params1.set<NonlinearVariableName>("variable") = var_name_v;
    _problem->addKernel("ConstantKernel", "ConstantKernel_" + var_name_v + "_" + Moose::stringify(cur_num), params1);
  }

  // second add mobile i
  for (int cur_num = 1; cur_num <= num_mobile_i; ++cur_num)
  {
    std::string var_name_i = name() + "0i" + Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("GMobile");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<std::vector<VariableName> > ("coupled_v_vars") = coupled_v_vars;
    params.set<std::vector<VariableName> > ("coupled_i_vars") = coupled_i_vars;
    params.set<UserObjectName>("user_object") = group_constant_name;
    params.set<unsigned int>("number_v") = number_v;
    params.set<unsigned int>("number_i") = number_i;
    params.set<unsigned int>("max_mobile_v") = num_mobile_v;
    params.set<unsigned int>("max_mobile_i") = num_mobile_i;
    _problem->addKernel("GMobile", "GMobile_" + var_name_i + "_" + Moose::stringify(cur_num), params);
    //printf("add GMobile: %s \n",var_name_i.c_str());

    // add pesudo kernel for L1 coefficient
    var_name_i = name() +"1i"+ Moose::stringify(cur_num);
    InputParameters params1 = _factory.getValidParams("ConstantKernel");
    params1.set<NonlinearVariableName>("variable") = var_name_i;
    _problem->addKernel("ConstantKernel", "ConstantKernel_" + var_name_i + "_" + Moose::stringify(cur_num), params1);
  }
}

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

#include "AddUserObjectDiffusion.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddUserObjectDiffusion, "add_kernel");

template<>
InputParameters validParams<AddUserObjectDiffusion>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::vector<unsigned int> >("mobile_v_size", "A vector of mobile species");
  params.addRequiredParam<std::vector<unsigned int> >("mobile_i_size", "A vector of mobile species");
  params.addRequiredParam<UserObjectName>("group_constant", "user object name");
  return params;
}

AddUserObjectDiffusion::AddUserObjectDiffusion(const InputParameters & params) :
    Action(params)
{
}

//only emission of point defects of same type are considered
void
AddUserObjectDiffusion::act()
{
  std::vector<unsigned int> v_size = getParam<std::vector<unsigned int> >("mobile_v_size");
  std::vector<unsigned int> i_size = getParam<std::vector<unsigned int> >("mobile_i_size");
  std::string group_constant_name = getParam<UserObjectName>("group_constant");
  std::vector<Real> vv, ii;

  const auto total_mobile_v = v_size.size();
  for (unsigned int i = 0; i < total_mobile_v; ++i)
    // account for sign; sink "+"
    vv.push_back(1.0);

  const auto total_mobile_i = i_size.size();
  for (unsigned int i = 0;i < total_mobile_i; ++i)
    //account for sign; sink "+"
    ii.push_back(1.0);

  for (unsigned int cur_num = 1; cur_num < total_mobile_v; ++cur_num)
  {
    // cur_num to vacancy, "+"
    Real coef = vv[cur_num-1];
    std::string var_name_v = name() + "v" + Moose::stringify(v_size[cur_num-1]);//add kernel for mobile v
    InputParameters params = _factory.getValidParams("UserObjectDiffusion");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    params.set<Real>("coeff") = coef; // loss "+"
    params.set<UserObjectName>("user_object") = group_constant_name;
    _problem->addKernel("UserObjectDiffusion", "Diffusion_" + var_name_v+Moose::stringify(cur_num), params);
    _console << "add UserObjectDiffusion disl: " << var_name_v << ", coef: " << coef << '\n';
  }

  for (unsigned int cur_num = 1; cur_num < total_mobile_i; ++cur_num)
  {
    Real coef = ii[cur_num-1];
    // add kernel for mobile i
    std::string var_name_i = name() + "i" + Moose::stringify(i_size[cur_num-1]);
    InputParameters params = _factory.getValidParams("UserObjectDiffusion");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    params.set<Real>("coeff") = coef; // loss "+"
    params.set<UserObjectName>("user_object") = group_constant_name;
    _problem->addKernel("UserObjectDiffusion", "Diffusion_" + var_name_i + Moose::stringify(cur_num), params);
    _console << "add UserObjectDiffusion disl: " << var_name_i << ", coef: " << coef << '\n';
  }
}

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

#include "GBase.h"
#include "Conversion.h"

template<>
InputParameters validParams<GBase>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("number_v", "Maximum vacancy cluster size");
  params.addRequiredParam<unsigned int>("number_i", "Maximum interstitial cluster size");
  params.addCoupledVar("coupled_v_vars", "coupled vacancy type variables");
  params.addCoupledVar("coupled_i_vars", "coupled intersitial type variables");
  params.addRequiredParam<unsigned int>("max_mobile_v", "A vector of mobile species");
  params.addRequiredParam<unsigned int>("max_mobile_i", "A vector of mobile species");
  params.addRequiredParam<UserObjectName>("user_object", "The name of user object providing interaction constants");
  return params;
}

GBase::GBase(const InputParameters & parameters)
     :Kernel(parameters),
     _number_v(getParam<unsigned int>("number_v")),
     _number_i(getParam<unsigned int>("number_i")),
     _max_mobile_v(getParam<unsigned int>("max_mobile_v")),
     _max_mobile_i(getParam<unsigned int>("max_mobile_i")),
     _gc(getUserObject<GGroup>("user_object"))
{
  NonlinearVariableName cur_var_name = getParam<NonlinearVariableName>("variable");
  _cur_size = getGroupNumber(cur_var_name);
  if (_cur_size > 0 && _cur_size < _max_mobile_v)
    mooseError("Please check the mobile vacancy cluster size range, should be immobile > mobile");
  else if (_cur_size < 0 && -_cur_size < _max_mobile_i)
    mooseError("Please check the mobile interstitial cluster size range, should be immobile>mobile");

  auto num_v_coupled = coupledComponents("coupled_v_vars");
  _no_v_vars.resize(num_v_coupled);
  _val_v_vars.resize(num_v_coupled);
  for (unsigned int i=0; i < num_v_coupled; ++i)
  {
    _no_v_vars[i] = coupled("coupled_v_vars", i);
    _val_v_vars[i] = &coupledValue("coupled_v_vars", i);
  }

  auto num_i_coupled = coupledComponents("coupled_i_vars");
  _no_i_vars.resize(num_i_coupled);
  _val_i_vars.resize(num_i_coupled);
  for (unsigned int i=0; i < num_i_coupled; ++i)
  {
    _no_i_vars[i] = coupled("coupled_i_vars", i);
    _val_i_vars[i] = &coupledValue("coupled_i_vars", i);
  }

#ifndef NDEBUG
  std::vector<VariableName> coupled_v_vars = getParam<std::vector<VariableName> >("coupled_v_vars");
  _console << "GBase: current variable => " << cur_var_name << std::endl;
  _console << "coupled with: " << std::endl;
  for (int i=0; i < num_v_coupled; ++i){
    _console << i << ":" << coupled_v_vars[i] << "  ";
  }
  _console << std::endl;
#endif
}

Real
GBase::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
  return 0.0;
}

int
GBase::getGroupNumber(std::string str)
{
  int len = str.length(), i = len;
  while (std::isdigit(str[i-1])) i--;

  int no = std::atoi((str.substr(i)).c_str());
  while (i>=0){
    i--;
    if (str[i]=='v')
      break;
    if (str[i]=='i')
    {
      no = -no;
      break;
    }
  }

  return no;
}

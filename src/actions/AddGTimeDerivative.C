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

#include "AddGTimeDerivative.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

template<>
InputParameters validParams<AddGTimeDerivative>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  return params;
}

AddGTimeDerivative::AddGTimeDerivative(const InputParameters & params) :
    Action(params)
{
}

void
AddGTimeDerivative::act()
{
  const auto number_v = getParam<unsigned int>("number_v");
  const auto number_i = getParam<unsigned int>("number_i");

  std::string var_name;
  for (unsigned int cur_num = 1; cur_num <= number_v; ++cur_num)
  {
    var_name = name() +"0v"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("TimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name;
    _problem->addKernel("TimeDerivative", "dt_0v_"+ var_name + Moose::stringify(cur_num), params);
    // printf("add TimeDerivative: %s\n",var_name_v.c_str());

    var_name = name() +"1v"+ Moose::stringify(cur_num);
    InputParameters params1 = _factory.getValidParams("TimeDerivative");
    params1.set<NonlinearVariableName>("variable") = var_name;
    _problem->addKernel("TimeDerivative", "dt_1v_" + var_name + Moose::stringify(cur_num), params1);
    // printf("add TimeDerivative: %s\n",var_name_v.c_str());
  }

  for (unsigned int cur_num = 1; cur_num <= number_i; cur_num++)
  {
    var_name = name() +"0i"+ Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("TimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name;
    _problem->addKernel("TimeDerivative", "dt_0i_"+ var_name + Moose::stringify(cur_num), params);
    // printf("add TimeDerivative: %s\n",var_name_i.c_str());

    var_name = name() +"1i"+ Moose::stringify(cur_num);
    InputParameters params1 = _factory.getValidParams("TimeDerivative");
    params1.set<NonlinearVariableName>("variable") = var_name;
    _problem->addKernel("TimeDerivative", "dt_1i_"+ var_name + Moose::stringify(cur_num), params1);
    // printf("add TimeDerivative: %s\n",var_name_i.c_str());
  }
}

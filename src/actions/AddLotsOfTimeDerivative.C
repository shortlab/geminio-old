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

#include "AddLotsOfTimeDerivative.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddLotsOfTimeDerivative, "add_kernel");

template<>
InputParameters validParams<AddLotsOfTimeDerivative>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  return params;
}

AddLotsOfTimeDerivative::AddLotsOfTimeDerivative(const InputParameters & params) :
    Action(params)
{
}

void
AddLotsOfTimeDerivative::act()
{
  const auto number_v = getParam<unsigned int>("number_v");
  const auto number_i = getParam<unsigned int>("number_i");

  for (unsigned int cur_num = 1; cur_num <= number_v; ++cur_num)
  {
    std::string var_name_v = name() + "v" + Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("TimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name_v;
    _problem->addKernel("TimeDerivative", "dt_v_" + var_name_v + Moose::stringify(cur_num), params);
  }

  for (unsigned int cur_num = 1; cur_num <= number_i; cur_num++)
  {
    std::string var_name_i = name() + "i" + Moose::stringify(cur_num);
    InputParameters params = _factory.getValidParams("TimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name_i;
    _problem->addKernel("TimeDerivative", "dt_i_" + var_name_i + Moose::stringify(cur_num), params);
  }
}

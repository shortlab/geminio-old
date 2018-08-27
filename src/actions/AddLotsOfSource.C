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

#include "AddLotsOfSource.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddLotsOfSource, "add_kernel");

template<>
InputParameters validParams<AddLotsOfSource>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<Real> >("source_v_size", "A vector of distribution of creation along ion range");
  params.addRequiredParam<std::vector<Real> >("source_i_size", "A vector of distribution of creation along ion range");
  params.addParam<std::string>("func_pre_name", "the sub-block name of [functions]");
  params.addParam<bool>("custom_input", false, "Use manually input coefficients");
  params.addParam<std::vector<Real> >("source_i", "production of i clusters for source_i_size species");
  params.addParam<std::vector<Real> >("source_v", "production of v clusters for 1 source_v_size species"); 
  return params;
}

AddLotsOfSource::AddLotsOfSource(const InputParameters & params) :
    Action(params)
{
  // if provided source_i or source_v (constant for each size), the function wouldn't be necessary. 
  // In other words, provide either function or source_i and source_v
  if (isParamValid("func_pre_name") && (isParamValid("source_i") || isParamValid("source_i")))
    paramError("func_pre_name", "Specify either 'func_pre_name' or 'source_i' and 'source_i'");
}

void
AddLotsOfSource::act()
{
  auto v_size = getParam<std::vector<Real> >("source_v_size");
  auto i_size = getParam<std::vector<Real> >("source_i_size");
  auto vv = getParam<std::vector<Real> >("source_v");
  auto ii = getParam<std::vector<Real> >("source_i");
  const auto custom = getParam<bool>("custom_input");
  const auto _total_v = getParam<unsigned int>("number_v");
  const auto _total_i = getParam<unsigned int>("number_i");
  std::string func_pre_name = (isParamValid("func_pre_name") ? getParam<std::string>("func_pre_name") : "");

  // ATTENTION: the emission of vacancy cluster emit an interstitial or interstitial cluster emit an vacancy is not considered
  for (unsigned int cur_num = 1; cur_num <= v_size.size(); cur_num++)
  {
    if (cur_num <= _total_v)
    {
      std::string var_name_v = name() +"v"+ Moose::stringify(v_size[cur_num-1]);
      std::string fun_name_v = func_pre_name + "v" + Moose::stringify(v_size[cur_num-1]);
      InputParameters params = _factory.getValidParams("BodyForce");
      params.set<NonlinearVariableName>("variable") = var_name_v;
      if (!custom){
        params.set<FunctionName>("function") = fun_name_v;
        params.set<Real>("value") = 1.0;
      }
      else if (vv.size() != v_size.size())
        mooseError("Vacancy source number should be the same with provided variables");
      else params.set<Real>("value") = vv[cur_num-1];//Should be the production term of current size, gain should be negative in the kernel
      _problem->addKernel("BodyForce", "bodyforce_" +  var_name_v + Moose::stringify(cur_num), params);
      _console << "add Source: " << var_name_v << '\n';
    }
  }
  for (unsigned int cur_num = 1; cur_num <= i_size.size(); cur_num++)
  {
    if (cur_num <= _total_i)
    {
      std::string var_name_i = name() +"i"+ Moose::stringify(i_size[cur_num-1]);
      std::string fun_name_i = func_pre_name + "i" + Moose::stringify(i_size[cur_num-1]);
      InputParameters params = _factory.getValidParams("BodyForce");
      params.set<NonlinearVariableName>("variable") = var_name_i;
      if (!custom){
        params.set<FunctionName>("function") = fun_name_i;
        params.set<Real>("value") = 1.0;
      }
      else if (ii.size() != i_size.size())
        mooseError("Intersitial source number should be the same with provided variable");
      else params.set<Real>("value") = ii[cur_num-1];// gain should be negative in the kernel
      _problem->addKernel("BodyForce", "bodyforce_" + var_name_i+ Moose::stringify(cur_num), params);
      _console << "add Source: " << var_name_i << '\n';
    }
  }
}

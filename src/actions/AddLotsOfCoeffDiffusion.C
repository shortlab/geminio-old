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

#include "AddLotsOfCoeffDiffusion.h"
#include "MaterialParameters.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Conversion.h"

registerMooseAction("GeminioApp", AddLotsOfCoeffDiffusion, "add_kernel");

template<>
InputParameters validParams<AddLotsOfCoeffDiffusion>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<unsigned int> >("mobile_v_size", "A vector of mobile species sizes");
  params.addRequiredParam<std::vector<unsigned int> >("mobile_i_size", "A vector of mobile species sizes");
  params.addParam<bool>("custom_input",false,"mannually input coefficients");
  params.addParam<std::vector<Real> >("diff_i", "diffusivity of i clusters from 1 to n_i");
  params.addParam<std::vector<Real> >("diff_v", "diffusivity of v clusters from 1 to n_v");
  params.addParam<Real>("temperature",800,"system temperature [K]");
  return params;
}

AddLotsOfCoeffDiffusion::AddLotsOfCoeffDiffusion(const InputParameters & params) :
    Action(params)
{
}

void
AddLotsOfCoeffDiffusion::act()
{
  //unsigned int number_v = getParam<unsigned int>("number_v");
  //unsigned int number_i = getParam<unsigned int>("number_i");
  const auto custom = getParam<bool>("custom_input");
  const auto _total_v = getParam<unsigned int>("number_v");
  const auto _total_i = getParam<unsigned int>("number_i");
  std::vector<unsigned int> v_size = getParam<std::vector<unsigned int> >("mobile_v_size");
  std::vector<unsigned int> i_size = getParam<std::vector<unsigned int> >("mobile_i_size");
  std::vector<Real> vv,ii;
  const auto temp = getParam<Real>("temperature");
  if (isParamValid("diff_v") && custom == true)
    vv = getParam<std::vector<Real> >("diff_v");
  else {
    for (unsigned int i = 0; i < v_size.size();i++){
      vv.push_back(diff(v_size[i],"V",temp));}
  }
  if (isParamValid("diff_i") && custom == true)
    ii = getParam<std::vector<Real> >("diff_i");
  else {
    for (unsigned int i = 0; i < i_size.size();i++){
      ii.push_back(diff(i_size[i],"I",temp));}
  }
  for (unsigned int cur_num = 1; cur_num <= v_size.size(); cur_num++)
  {
    if (cur_num <= _total_v){
      std::string var_name_v = name() +"v"+ Moose::stringify(v_size[cur_num-1]);
      InputParameters params0 = _factory.getValidParams("CoeffDiffusion");
      params0.set<NonlinearVariableName>("variable") = var_name_v;
      params0.set<Real>("diffusivity") =  vv[cur_num-1];
      _problem->addKernel("CoeffDiffusion", "diffusion_" + var_name_v + "_" + Moose::stringify(cur_num), params0);
    }
  }
  for (unsigned int cur_num = 1; cur_num <= i_size.size(); cur_num++)
  {
    if (cur_num <= _total_i){
      std::string var_name_i = name() + "i" + Moose::stringify(i_size[cur_num-1]);
      InputParameters params0 = _factory.getValidParams("CoeffDiffusion");
      params0.set<NonlinearVariableName>("variable") = var_name_i;
      params0.set<Real>("diffusivity") =  vv[cur_num-1];
      _problem->addKernel("CoeffDiffusion", "diffusion_" + var_name_i + '_' + Moose::stringify(cur_num), params0);
    }
  }
}

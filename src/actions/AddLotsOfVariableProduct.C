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

#include "AddLotsOfVariableProduct.h"
#include "MaterialParameters.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "VariableProduct.h"

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
InputParameters validParams<AddLotsOfVariableProduct>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<unsigned int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<unsigned int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<unsigned int> >("mobile_v_size", "A vector of mobile species sizes");
  params.addRequiredParam<std::vector<unsigned int> >("mobile_i_size", "A vector of mobile species sizes");
  params.addParam<Real>("temperature",600.0,"Temperature");
  params.addParam<bool>("custom_input",false,"true: use manually input data");
  params.addParam<std::vector<Real> >("absorb_vv", "absorption coefficient between v and v; row of 1 to n_v and column of 1 to n_v (n_v x n_v)");
  params.addParam<std::vector<Real> >("absorb_ii", "absorption coefficient between v and v; row of 1 to n_i and column of 1 to n_i (n_i x n_i)");
  params.addParam<std::vector<Real> >("absorb_vi", "absorption coefficient between i and v; row of 1 to n_v and column of 1 to n_i (n_v x n_i)");
  return params;
}

AddLotsOfVariableProduct::AddLotsOfVariableProduct(const InputParameters & params) :
    AddVariableAction(params)
{
}

void
AddLotsOfVariableProduct::act()
{
  const auto number_v = getParam<unsigned int>("number_v");
  const auto number_i = getParam<unsigned int>("number_i");
  std::vector<unsigned int> v_size = getParam<std::vector<unsigned int> >("mobile_v_size");
  std::vector<unsigned int> i_size = getParam<std::vector<unsigned int> >("mobile_i_size");
  const auto temp = getParam<Real>("temperature");
  const auto custom = getParam<bool>("custom_input");
  std::vector<Real> vv, ii, vi;

  // used to mark mobility, 1: mobile; 0: immobile
  int tagi = 0, tagj = 0;

  if (isParamValid("absorb_vv") && custom == true)
    vv = getParam<std::vector<Real> >("absorb_vv");
  else {
    for (unsigned int i = 1; i <= number_v; ++i)
    {
      if (std::find(v_size.begin(),v_size.end(),i) != v_size.end())
        tagi = 1;
      else
        tagi = 0;

      for (int j=1;j<=number_v;j++)
      {
        if (std::find(v_size.begin(),v_size.end(),j) != v_size.end())//confirm it as mobile
          tagj = 1;
        else
          tagj = 0;

        // TODO: need to add mobile information to the function
        vv.push_back(absorb(i, j, "V", "V", temp, tagi, tagj));
        // printf("absorb v %d with v %d: %lf\n",i,j,vv.back());
      }
    }
  }

  if (isParamValid("absorb_ii") && custom == true)
    ii = getParam<std::vector<Real> >("absorb_ii");
  else {
    for (unsigned int i = 1; i <= number_i; ++i)
    {
      if (std::find(i_size.begin(),i_size.end(),i) != i_size.end())
        tagi = 1; //mobile
      else
        tagi = 0; //immobile
      for (unsigned int j = 1; j <= number_i; ++j)
      {
        if (std::find(i_size.begin(),i_size.end(),j) != i_size.end())//confirm it as mobile
          tagj = 1; // mobile
        else
          tagj = 0; // immobile
        ii.push_back(absorb(i, j, "I", "I", temp, tagi, tagj));
        // printf("absorb i %d with i %d: %lf\n",i,j,ii.back());
      }
    }
  }

  if (isParamValid("absorb_vi") && custom == true)
    vi = getParam<std::vector<Real> >("absorb_vi");
  else {
    for (unsigned int i = 1; i <= number_v; ++i){
      if (std::find(v_size.begin(),v_size.end(),i) != v_size.end())
        tagi = 1; // mobile
      else
        tagi = 0; // immobile

      for (unsigned int j = 1; j <= number_i; ++j)
      {
        // confirm it as mobile
        if (std::find(i_size.begin(),i_size.end(),j) != i_size.end())
          tagj = 1; // mobile
        else
          tagj = 0; // immobile
        vi.push_back(absorb(i, j, "V", "I", temp, tagi, tagj));
        // printf("absorb v %d with i %d: %lf\n",i,j,vi.back());
      }
    }
  }

  // First comes with the kernels for vacancy clusters
  for (int cur_num = 1; cur_num <= number_v; cur_num++)
  {
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);

    // vi reaction down(-)
    for (int cur_num2 = 1; cur_num2 <= number_i; cur_num2++)
    {
      Real coef = 0.0;

      // confirm it as mobile
      if (std::find(i_size.begin(),i_size.end(),cur_num2) != i_size.end() ||
          std::find(v_size.begin(),v_size.end(),cur_num) != v_size.end() )
      {
        std::string var_name2;
        var_name2 = name() +"i"+ Moose::stringify(cur_num2);
        InputParameters params1 = _factory.getValidParams("VariableProduct");
        params1.set<NonlinearVariableName>("variable") = var_name_v;
        params1.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
        coef =  vi[(cur_num-1)*number_i+cur_num2-1];
        params1.set<Real>("coeff") = coef;
        _problem->addKernel("VariableProduct", "VarProd_" + var_name_v + "_vim_" + Moose::stringify(cur_num2), params1);
        //printf("add VariableProduct: %s (%s and %s), coef: %lf\n",var_name_v.c_str(),var_name_v.c_str(),var_name2.c_str(),coef);
      }
    }

    // vv reaction down(-)
    // ATTENTION: largest number_v-cur_num to ensure conservation (no outlier reaction)
    for (int cur_num2 = 1; cur_num2 <= number_v-cur_num; cur_num2++)
    {
      Real coef = 0.0;

      // confirm it as mobile
      if (std::find(v_size.begin(),v_size.end(),cur_num2) != v_size.end() ||
          std::find(v_size.begin(),v_size.end(),cur_num) != v_size.end())
      {
        std::string var_name2;
        var_name2 = name() +"v"+ Moose::stringify(cur_num2);
        InputParameters params1 = _factory.getValidParams("VariableProduct");
        params1.set<NonlinearVariableName>("variable") = var_name_v;
        if (cur_num2 == cur_num){//if the same species, then no need to add second variable, the coefficient is multiplied by 2
          coef = 2*vv[(cur_num-1)*number_v+cur_num2-1];
          params1.set<Real>("coeff") = coef;
        }
        else{
          params1.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
          coef = vv[(cur_num-1)*number_v+cur_num2-1];
          params1.set<Real>("coeff") = coef;
        }
        _problem->addKernel("VariableProduct", "VarProd_" + var_name_v+ "_vvm_" + Moose::stringify(cur_num2), params1);
        //printf("add VariableProduct: %s (%s and %s), coef: %lf\n",var_name_v.c_str(),var_name_v.c_str(),var_name2.c_str(),coef);
      }
    }
  }

  // int pairs_v = (number_v + 1) / 2;

  for (unsigned int cur_num = 1; cur_num <= number_v; cur_num++)
  {
    std::string var_name_v = name() +"v"+ Moose::stringify(cur_num);
    //for (int j=1; j<= pairs_v; j++)

    // vv reaction up(+)
    for (unsigned int j = 1; j <= cur_num / 2; ++j)
    {
      Real coef = 0.0;

      // at least one is mobile
      if (std::find(v_size.begin(),v_size.end(),j) != v_size.end() ||
          std::find(v_size.begin(),v_size.end(),cur_num-j) != v_size.end())
      {
        std::string var_name2 = name() +"v"+ Moose::stringify(j);
        std::string var_name3 = name() +"v"+ Moose::stringify(cur_num-j);
        InputParameters params = _factory.getValidParams("VariableProduct");
        params.set<NonlinearVariableName>("variable") = var_name_v;
        params.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
        params.set<std::vector<VariableName> >("coupled_vars").push_back(var_name3);
        coef = -vv[(j-1) * number_v + cur_num - j - 1];

        //gain should be negative in the kernel
        params.set<Real>("coeff") = coef;

        _problem->addKernel("VariableProduct", "VarProd_" + var_name_v+ "_vvp_" + Moose::stringify(j), params);
        //printf("add VariableProduct: %s (%s and %s), coef: %lf\n",var_name_v.c_str(),var_name2.c_str(),var_name3.c_str(),coef);
      }
    }

    // vi reaction up(+)
    for (unsigned int j = cur_num + 1; j<= number_v; ++j)
    {
      Real coef = 0.0;

      // at least one is mobile
      if (std::find(v_size.begin(),v_size.end(),j) != v_size.end() ||
          std::find(i_size.begin(),i_size.end(),j-cur_num) != i_size.end())
      {
        // guarantee the maximum size of interstitial clusters satisfy
        if (j-cur_num <= number_i)
        {
          std::string var_name2 = name() + "v" + Moose::stringify(j);
          std::string var_name3 = name() + "i" + Moose::stringify(j-cur_num);
          InputParameters params = _factory.getValidParams("VariableProduct");
          params.set<NonlinearVariableName>("variable") = var_name_v;
          params.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
          params.set<std::vector<VariableName> >("coupled_vars").push_back(var_name3);
          coef = -vi[(j-1) * number_i + j - cur_num -1];

          //gain should be negative in the kernel, factor 2 account for mutual position
          params.set<Real>("coeff") = coef;

          _problem->addKernel("VariableProduct", "VarProd_" + var_name_v+ "_vip_" + Moose::stringify(j), params);
          //printf("add VariableProduct: %s (%s and %s), coef: %lf\n",var_name_v.c_str(),var_name2.c_str(),var_name3.c_str(),coef);
        }
      }
    }
  }

  // Second comes with the kernels for interstitial clusters
  for (unsigned int cur_num = 1; cur_num <= number_i; ++cur_num)
  {
    std::string var_name_i = name() + "i" + Moose::stringify(cur_num);

    // ii reaction down(-)
    // ATTENTION: largest number_v-cur_num to ensure conservation (no outlier reaction)
    for (unsigned int cur_num2 = 1; cur_num2 <= number_i - cur_num; ++cur_num2)
    {
      Real coef = 0.0;

      // confirm it as mobile
      if (std::find(i_size.begin(),i_size.end(),cur_num2) != i_size.end() ||
          std::find(i_size.begin(),i_size.end(),cur_num) != i_size.end() )
      {
        std::string var_name2;
        var_name2 = name() +"i"+ Moose::stringify(cur_num2);
        InputParameters params1 = _factory.getValidParams("VariableProduct");
        params1.set<NonlinearVariableName> ("variable") = var_name_i;
        if (cur_num == cur_num2)
        {
          //same species, factor 2 to account for mutual position
          coef = 2*ii[(cur_num-1)*number_i+cur_num2-1];
          params1.set<Real>("coeff") = coef;
        }
        else
        {
          params1.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
          coef = ii[(cur_num-1)*number_i+cur_num2-1];
          params1.set<Real>("coeff") = coef;
        }

        _problem->addKernel("VariableProduct", "VarProd_" + var_name_i + "_iim_" + Moose::stringify(cur_num2), params1);
        //printf("add VariableProduct: %s (%s and %s), coef: %lf\n",var_name_i.c_str(),var_name_i.c_str(),var_name2.c_str(),coef);
      }
    }

    // iv reaction down(-)
    for (unsigned int cur_num2 = 1; cur_num2 <= number_v; ++cur_num2)
    {
      Real coef = 0.0;

      // confirm it as mobile
      if (std::find(v_size.begin(),v_size.end(),cur_num2) != v_size.end() ||
          std::find(i_size.begin(),i_size.end(),cur_num) != i_size.end())
      {
        std::string var_name2;
        var_name2 = name() + "v" + Moose::stringify(cur_num2);
        InputParameters params1 = _factory.getValidParams("VariableProduct");
        params1.set<NonlinearVariableName> ("variable") = var_name_i;
        params1.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
        coef = vi[(cur_num2-1)*number_i+cur_num-1];
        params1.set<Real>("coeff") = coef;
        _problem->addKernel("VariableProduct", "VarProd_" + var_name_i + "_ivm_" + Moose::stringify(cur_num2), params1);
        //printf("add VariableProduct: %s (%s and %s), coef: %lf\n",var_name_i.c_str(),var_name_i.c_str(),var_name2.c_str(),coef);
      }
    }
  }

  // int  pairs_i = (number_i+1)/2;
  for (unsigned int cur_num = 1; cur_num <= number_i; ++cur_num)
  {
    std::string var_name_i = name() + "i" + Moose::stringify(cur_num);
    //for (int j=1; j<= pairs_i; j++)

    // ii reaction up(+)
    for (unsigned int j = 1; j <= cur_num / 2; ++j)
    {
      Real coef = 0.0;
      // at least one is mobile
      if (std::find(i_size.begin(),i_size.end(),j) != i_size.end() ||
          std::find(i_size.begin(),i_size.end(),cur_num-j) != i_size.end())
      {
          std::string var_name2 = name() +"i"+ Moose::stringify(j);
          std::string var_name3 = name() +"i"+ Moose::stringify(cur_num-j);
          InputParameters params = _factory.getValidParams("VariableProduct");
          params.set<NonlinearVariableName>("variable") = var_name_i;
          params.set<std::vector<VariableName> > ("coupled_vars").push_back(var_name2);
          params.set<std::vector<VariableName> >("coupled_vars").push_back(var_name3);
          coef = (-1.0)*ii[(j-1)*number_i+cur_num-j-1];
          // gain should be negative in the kernel
          params.set<Real>("coeff") = coef;
          _problem->addKernel("VariableProduct", "VarProd_" + var_name_i+ "_iip_" + Moose::stringify(j), params);
          //printf("add VariableProduct: %s (%s and %s), coef: %lf\n",var_name_i.c_str(),var_name2.c_str(),var_name3.c_str(),coef);
      }
    }

    //iv reaction up(+)
    for (unsigned int j = cur_num + 1; j <= number_i; ++j)
    {
      Real coef = 0.0;
      //at least one is mobile
      if (std::find(i_size.begin(),i_size.end(),j) != i_size.end() ||
          std::find(v_size.begin(),v_size.end(),j-cur_num) != v_size.end())
      {
        // guarantee the maximum size of interstitial clusters satisfy
        if (j-cur_num <= number_v)
        {
          std::string var_name2 = name() +"i"+ Moose::stringify(j);
          std::string var_name3 = name() +"v"+ Moose::stringify(j-cur_num);
          InputParameters params = _factory.getValidParams("VariableProduct");
          params.set<NonlinearVariableName>("variable") = var_name_i;
          params.set<std::vector<VariableName> >("coupled_vars").push_back(var_name2);
          params.set<std::vector<VariableName> >("coupled_vars").push_back(var_name3);
          coef = (-1.0)*vi[(j-cur_num-1)*number_i+j-1];
          // gain should be negative in the kernel
          params.set<Real>("coeff") = coef;
          _problem->addKernel("VariableProduct", "VarProd_" + var_name_i + "_ivp_" + Moose::stringify(j), params);
          //printf("add VariableProduct: %s (%s and %s), coef: %lf\n",var_name_i.c_str(),var_name2.c_str(),var_name3.c_str(),coef);
        }
      }
    }
  }
}

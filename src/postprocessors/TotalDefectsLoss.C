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

///////////////////////// Check to total number of defects in terms of point defect /////////////

#include "TotalDefectLoss.h"
#include "MooseMesh.h"
#include "SubProblem.h"
#include "Conversion.h"

// libMesh
#include "libmesh/node.h"

registerMooseObject("GeminioApp", TotalDefectLoss);

template<>
InputParameters validParams<TotalDefectLoss>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<unsigned int>("nodeid", "The ID of the node where we monitor");
  params.addRequiredParam<std::string>("var_prefix_v", "The prefix string of variables (in front of number)");
  params.addRequiredParam<std::string>("var_prefix_i", "The prefix string of variables (in front of number)");
  params.addRequiredParam<int>("number_v", "The number of vacancy variables to add");
  params.addRequiredParam<int>("number_i", "The number of interstitial variables to add");
  params.addRequiredParam<std::vector<int> >("mobile_v_size", "A vector of mobile species sizes");
  params.addRequiredParam<std::vector<int> >("mobile_i_size", "A vector of mobile species sizes");
  params.addRequiredParam<UserObjectName>("material", "User object to retrieve interaction coefficient");
  params.addRequiredParam<Real>("temperature", "system temperature");
  return params;
}

TotalDefectLoss::TotalDefectLoss(const InputParameters & parameters) 
  : GeneralPostprocessor(parameters),
    _mesh(_subproblem.mesh()),
    _var_prefix_v(getParam<std::string>("var_prefix_v")),
    _var_prefix_i(getParam<std::string>("var_prefix_i")),
    _mobile_v(getParam<std::vector<int> >("mobile_v_size")),
    _mobile_i(getParam<std::vector<int> >("mobile_i_size")),
    _v_max(getParam<int>("number_v")),
    _i_max(getParam<int>("number_i")),
    _total_loss(declareRestartableData<Real>("TotalLoss", 0.0)),
    _node_ptr(_mesh.getMesh().query_node_ptr(getParam<unsigned int>("nodeid"))),
    _material(&getUserObject<MaterialConstants>("material")),
    _T(getParam<Real>("temperature"))
{
  // This class only works with ReplicatedMesh, since it relies on a
  // specific node numbering that we can't guarantee with DistributedMesh
  _mesh.errorIfDistributedMesh("TotalDefectLoss");

  if (_node_ptr == nullptr)
    paramError("nodeid", "Node #", getParam<unsigned int>("nodeid"), " not found in the mesh!");
}

Real
TotalDefectLoss::getValue()
{
  Real coef = 0.0;
  Real value = 0.0;
  std::string var_name_v, var_name_i;
  int loss_size, tagi, tagj;
  int v_size = _mobile_v.size();
  int i_size = _mobile_i.size();
  
  tagj = 1;
  for (int vv=1;vv<=_v_max;vv++)
  {
    tagi = 0;

    // at least one is mobile
    if (std::find(_mobile_v.begin(), _mobile_v.end(), vv) != _mobile_v.end())
      tagi = 1;

    for (int ii = 0; ii < i_size; ++ii)
    {
      value = 0.0;
      var_name_v = _var_prefix_v + Moose::stringify(vv);
      var_name_i = _var_prefix_i + Moose::stringify(_mobile_i[ii]);

      // get the smaller one
      loss_size = vv > _mobile_i[ii] ? _mobile_i[ii] : vv; 
      
      // absorption between i and j; TODO:need change
      coef = _material->absorb(vv, _mobile_i[ii], MaterialParameters::Species::V, MaterialParameters::Species::I, _T, tagi, tagj);
      
      if (_node_ptr->processor_id() == processor_id())
      {
        // defect size 
        const auto size_v = _subproblem.getVariable(_tid, var_name_v).getNodalValue(*_node_ptr);
        const auto size_i = _subproblem.getVariable(_tid, var_name_i).getNodalValue(*_node_ptr);
        value = 2.0 *_fe_problem.dt() * coef * size_v * size_i * loss_size; 
      } 

      // broadcast one value to other processors with previous zeros
      gatherSum(value);
      _total_loss += value;
    } 
  }
  
  tagi = 1;
  for (int ii = 1; ii <= _i_max; ++ii)
  {
    tagj = 0; 
    // at least one is mobile
    if (std::find(_mobile_i.begin(), _mobile_i.end(), ii) != _mobile_i.end())
      tagj = 1;
    
    for (int vv = 0; vv < v_size; ++vv)
    {
      value = 0.0;
      var_name_v = _var_prefix_v + Moose::stringify(_mobile_v[vv]);
      var_name_i = _var_prefix_i + Moose::stringify(ii);
      
      // get the smaller one
      loss_size = ii > _mobile_v[vv] ? _mobile_i[vv] : ii; 

      // absorption between i and j; TODO:need change
      coef = _material->absorb(ii, _mobile_v[vv], MaterialParameters::Species::I, MaterialParameters::Species::V, _T, tagi, tagj);
      if (_node_ptr->processor_id() == processor_id())
      {
        // defect size 
        const auto size_v = _subproblem.getVariable(_tid, var_name_v).getNodalValue(*_node_ptr);
        const auto size_i = _subproblem.getVariable(_tid, var_name_i).getNodalValue(*_node_ptr);
        value = 2.0* _fe_problem.dt() * coef * size_v * size_i * loss_size;
      }

      // broadcast one value to other processors with previous zeros
      gatherSum(value);
      _total_loss += value;
    } 
  }

  // subtract the double account of mobile-mobile reaction
  tagi = 1; tagj = 1;
  for (int ii = 0; ii < i_size; ++ii)
  {
    for (int vv = 0; vv < v_size; ++vv)
    {
      value = 0.0;
      var_name_v = _var_prefix_v + Moose::stringify(_mobile_v[vv]);
      var_name_i = _var_prefix_i + Moose::stringify(_mobile_i[ii]);
      
      // get the smaller one
      loss_size = _mobile_i[ii] > _mobile_v[vv] ? _mobile_i[vv] : _mobile_i[ii]; 
      
      // absorption between i and j; TODO:need change
      coef = _material->absorb(_mobile_i[ii], _mobile_v[vv], MaterialParameters::Species::I, MaterialParameters::Species::V, _T, tagi, tagj);
      if (_node_ptr->processor_id() == processor_id())
      {
        // defect size 
        const auto size_v = _subproblem.getVariable(_tid, var_name_v).getNodalValue(*_node_ptr);
        const auto size_i = _subproblem.getVariable(_tid, var_name_i).getNodalValue(*_node_ptr);
        value = -2.0* _fe_problem.dt()*coef * size_v * size_i * loss_size;
      }

      // broadcast one value to other processors with previous zeros
      gatherSum(value);
      _total_loss += value;
    } 
  }

  return _total_loss;
}


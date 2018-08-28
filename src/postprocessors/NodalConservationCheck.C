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

#include "NodalConservationCheck.h"
#include "MooseMesh.h"
#include "SubProblem.h"
#include "Conversion.h"

// libMesh
#include "libmesh/node.h"

registerMooseObject("GeminioApp", NodalConservationCheck);

template<>
InputParameters validParams<NodalConservationCheck>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<unsigned int>("nodeid", "The ID of the node where we monitor");
  params.addRequiredParam<std::string>("var_prefix", "The prefix string of variables (in front of number)");
  params.addRequiredParam<std::vector<int> >("size_range", "The number of vacancy variables to solve");
  params.addParam<Real>("scale_factor", 1, "A scale factor to be applied to the variable");
  return params;
}

NodalConservationCheck::NodalConservationCheck(const InputParameters & parameters) :
    GeneralPostprocessor(parameters),
    _mesh(_subproblem.mesh()),
    _var_prefix(getParam<std::string>("var_prefix")),
    _size_range(getParam<std::vector<int> >("size_range")),
    _node_ptr(_mesh.getMesh().query_node_ptr(getParam<unsigned int>("nodeid"))),
    _scale_factor(getParam<Real>("scale_factor"))
{
  // This class only works with ReplicatedMesh, since it relies on a
  // specific node numbering that we can't guarantee with DistributedMesh
  _mesh.errorIfDistributedMesh("NodalConservationCheck");

  if (_node_ptr == nullptr)
    paramError("nodeid", "Node #", getParam<unsigned int>("nodeid"), " not found in the mesh!");

  // neither provided or both provided is wrong
  if (_size_range.size() != 2 || _size_range[1] < _size_range[0])
    mooseError("Defect size range is not provided correctly, double check!");
}

Real
NodalConservationCheck::getValue()
{
  Real total = 0.0;
  int max_num = _size_range[1];

  Real value;
  for (int i = _size_range[0]; i <= max_num; ++i)
  {
    value = 0.0;
    std::string _var_name = _var_prefix + Moose::stringify(i);
    if (_node_ptr->processor_id() == processor_id())
    {
      const Real size = _subproblem.getVariable(_tid, _var_name).getNodalValue(*_node_ptr);
      value = size * i;
    }

    // broadcast one value to other processors with previous zeros
    gatherSum(value);

    total += value;
  }

  return _scale_factor * total;
}

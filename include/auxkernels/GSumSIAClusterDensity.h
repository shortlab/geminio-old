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

#ifndef GSUMSIACLUSTERDENSITY_H
#define GSUMSIACLUSTERDENSITY_H

#include "AuxKernel.h"
#include "GGroup.h"

//Forward Declarations
class GSumSIAClusterDensity;

template<>
InputParameters validParams<GSumSIAClusterDensity>();

/**
 * Calculate total density of SIA cluster in range [lower_bounds,upper_bound]
 * from variable values with grouping method
 */
class GSumSIAClusterDensity : public AuxKernel
{
public:
  GSumSIAClusterDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  const GGroup & _gc;
  Real _scale_factor;
  unsigned int _lower_bound;
  unsigned int _upper_bound;
  std::vector<unsigned int> _no_vars;
  std::vector<const VariableValue *> _val_vars;
};

#endif // GSUMSIACLUSTERDENSITY_H

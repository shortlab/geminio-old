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

#ifndef CLUSTERDENSITY_H
#define CLUSTERDENSITY_H

#include "AuxKernel.h"

//Forward Declarations
class ClusterDensity;

template<>
InputParameters validParams<ClusterDensity>();

/**
 * calculate void swelling from variable values with discrete cluster dynamics method
 */
class ClusterDensity : public AuxKernel
{
public:
  ClusterDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  Real _scaling_factor;
  std::vector<unsigned int> _no_vars;
  std::vector<const VariableValue *> _val_vars;
};

#endif // CLUSTERDENSITY_H

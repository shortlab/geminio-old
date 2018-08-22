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

#ifndef GVOIDSWELLING_H
#define GVOIDSWELLING_H

#include "AuxKernel.h"
#include "GGroup.h"

//Forward Declarations
class GVoidSwelling;

template<>
InputParameters validParams<GVoidSwelling>();

/**
 * Calculate void swelling from variable vaules with grouping method
 */
class GVoidSwelling : public AuxKernel
{
public:
  GVoidSwelling(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  const GGroup & _gc;
  std::vector<unsigned int> _no_v_vars;
  std::vector<const VariableValue *> _val_v_vars;
};

#endif // GVOIDSWELLING_H

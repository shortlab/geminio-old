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

#ifndef DEFECTSOURCE_H
#define DEFECTSOURCE_H

#include "Kernel.h"

//Forward Declarations
class DefectSource;

template<>
InputParameters validParams<DefectSource>();

class DefectSource : public Kernel
{
public:
  DefectSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  Real _coef;
  const VariableValue & _primary_source;
  const VariableValue & _secondary_source;
  const MaterialProperty<Real> & _intracascade_survival;
};

#endif // DEFECTSOURCE_H

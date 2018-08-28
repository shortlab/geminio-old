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

#ifndef GIMMOBILEL0_H
#define GIMMOBILEL0_H

#include "GBase.h"
#include "GGroup.h"

//Forward Declarations
class GImmobileL0;

template<>
InputParameters validParams<GImmobileL0>();

class GImmobileL0 : public GBase
{
public:
  GImmobileL0(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};

#endif // GIMMOBILEL0_H

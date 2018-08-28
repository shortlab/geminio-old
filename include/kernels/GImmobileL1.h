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

#ifndef GIMMOBILEL1_H
#define GIMMOBILEL1_H

#include "GBase.h"
#include "GGroup.h"

//Forward Declarations
class GImmobileL1;

template<>
InputParameters validParams<GImmobileL1>();

class GImmobileL1 : public GBase
{
public:
  GImmobileL1(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};

#endif // GIMMOBILEL1_H

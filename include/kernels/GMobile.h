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

#ifndef GMOBILE_H
#define GMOBILE_H

#include "GBase.h"
#include "GGroup.h"

//Forward Declarations
class GMobile;

template<>
InputParameters validParams<GMobile>();

class GMobile : public GBase
{
public:
  GMobile(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  double getConcBySize(int i);

protected:
  int _max_v;
  int _max_i;
};

#endif // GMOBILE_H

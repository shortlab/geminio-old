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

#ifndef ADDGTIMEDERIVATIVE_H
#define ADDGTIMEDERIVATIVE_H

#include "Action.h"

class AddGTimeDerivative;

template<>
InputParameters validParams<AddGTimeDerivative>();


class AddGTimeDerivative : public Action
{
public:
  AddGTimeDerivative(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDGTIMEDERIVATIVE_H

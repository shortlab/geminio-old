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

#ifndef ADDLOTSOFTIMEDERIVATIVE_H
#define ADDLOTSOFTIMEDERIVATIVE_H

#include "Action.h"

class AddLotsOfTimeDerivative;

template<>
InputParameters validParams<AddLotsOfTimeDerivative>();


class AddLotsOfTimeDerivative : public Action
{
public:
  AddLotsOfTimeDerivative(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDLOTSOFTIMEDERIVATIVE_H

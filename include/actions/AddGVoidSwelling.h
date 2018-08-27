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

#ifndef ADDGVOIDSWELLING_H
#define ADDGVOIDSWELLING_H

#include "Action.h"

class AddGVoidSwelling;

template<>
InputParameters validParams<AddGVoidSwelling>();


class AddGVoidSwelling : public Action
{
public:
  AddGVoidSwelling(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDGVOIDSWELLING_H

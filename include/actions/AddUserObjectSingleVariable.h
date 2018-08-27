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

#ifndef ADDUSEROBJECTSINGLEVARIABLE_H
#define ADDUSEROBJECTSINGLEVARIABLE_H

#include "Action.h"

class AddUserObjectSingleVariable;

template<>
InputParameters validParams<AddUserObjectSingleVariable>();


class AddUserObjectSingleVariable : public Action
{
public:
  AddUserObjectSingleVariable(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDUSEROBJECTSINGLEVARIABLE_H

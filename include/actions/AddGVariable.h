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

#ifndef ADDGVARIABLE_H
#define ADDGVARIABLE_H

#include "GeminioAddVariableAction.h"

class AddGVariable;

template<>
InputParameters validParams<AddGVariable>();


class AddGVariable : public GeminioAddVariableAction
{
public:
  AddGVariable(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDLOTSOFVARIABLEACTION_H

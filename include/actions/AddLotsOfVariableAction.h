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

#ifndef ADDLOTSOFVARIABLEACTION_H
#define ADDLOTSOFVARIABLEACTION_H

#include "GeminioAddVariableAction.h"

class AddLotsOfVariableAction;

template<>
InputParameters validParams<AddLotsOfVariableAction>();


class AddLotsOfVariableAction : public GeminioAddVariableAction
{
public:
  AddLotsOfVariableAction(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDLOTSOFVARIABLEACTION_H

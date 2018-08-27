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

#ifndef GEMINIOADDVARIABLEACTION_H
#define GEMINIOADDVARIABLEACTION_H

#include "AddVariableAction.h"
#include "MooseEnum.h"

class GeminioAddVariableAction;

template<>
InputParameters validParams<GeminioAddVariableAction>();


class GeminioAddVariableAction : public AddVariableAction
{
public:
  GeminioAddVariableAction(const  InputParameters & parameters);

protected:
  void addConstantIC(const std::string & var_name, Real initial);
  
  const MooseEnum _bc_type;
  std::string _bc_name;
  const Real _boundary_value;
  const unsigned int _number_i;
  const unsigned int _number_v;
};

#endif // GEMINIOADDVARIABLEACTION_H

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

#ifndef ADDUSEROBJECTVARIABLEPRODUCT_H
#define ADDUSEROBJECTVARIABLEPRODUCT_H

#include "Action.h"

class AddUserObjectVariableProduct;

template<>
InputParameters validParams<AddUserObjectVariableProduct>();


class AddUserObjectVariableProduct : public Action
{
public:
  AddUserObjectVariableProduct(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDUSEROBJECTVARIABLEPRODUCT_H

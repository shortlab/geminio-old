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

#ifndef ADDLOTSOFVARIABLEPRODUCT_H
#define ADDLOTSOFVARIABLEPRODUCT_H

#include "Action.h"

class AddLotsOfVariableProduct;

template<>
InputParameters validParams<AddLotsOfVariableProduct>();


class AddLotsOfVariableProduct : public Action
{
public:
  AddLotsOfVariableProduct(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDLOTSOFVARIABLEPRODUCT_H

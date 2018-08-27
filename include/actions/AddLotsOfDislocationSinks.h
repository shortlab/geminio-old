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

#ifndef ADDLOTSOFDISLOCATIONSINKS_H
#define ADDLOTSOFDISLOCATIONSINKS_H

#include "Action.h"

class AddLotsOfDislocationSinks;

template<>
InputParameters validParams<AddLotsOfDislocationSinks>();


class AddLotsOfDislocationSinks : public Action
{
public:
  AddLotsOfDislocationSinks(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDLOTSOFDISLOCATIONSINKS_H

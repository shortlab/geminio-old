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

#ifndef ADDGSUMSIACLUSTERDENSITY_H
#define ADDGSUMSIACLUSTERDENSITY_H

#include "Action.h"

class AddGSumSIAClusterDensity;

template<>
InputParameters validParams<AddGSumSIAClusterDensity>();


class AddGSumSIAClusterDensity : public Action
{
public:
  AddGSumSIAClusterDensity(const  InputParameters & parameters);

  virtual void act();
};

#endif // ADDGSUMSIACLUSTERDENSITY_H

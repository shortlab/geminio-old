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

#ifndef GBASE_H
#define GBASE_H

#include "Kernel.h"
#include "GGroup.h"

//Forward Declarations
class GBase;


template<>
InputParameters validParams<GBase>();

class GBase : public Kernel
{
public:
  GBase(const InputParameters & parameters);

protected:
  virtual Real computeQpOffDiagJacobian(unsigned int /*jvar*/);
  int getGroupNumber(std::string);

protected:
  unsigned int _number_v;
  unsigned int _number_i;
  int _max_mobile_v;
  int _max_mobile_i;
  const GGroup & _gc;
  std::vector<unsigned int> _no_v_vars;
  std::vector<const VariableValue *> _val_v_vars;
  std::vector<unsigned int> _no_i_vars;
  std::vector<const VariableValue *> _val_i_vars;
  int _cur_size;
};

#endif // GBASE_H

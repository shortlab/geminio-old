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

#ifndef DEFECTRECOMBINATION_H
#define DEFECTRECOMBINATION_H

#include "Kernel.h"

//Forward Declarations
class DefectRecombination;

template<>
InputParameters validParams<DefectRecombination>();

class DefectRecombination : public Kernel
{
public:
  DefectRecombination(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableValue & _other_defect;
  unsigned int _other_defect_var;
  const VariableValue & _recombination_rate;
};

#endif // DEFECTRECOMBINATION_H

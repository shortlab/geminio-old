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

#include "DefectRecombination.h"

registerMooseObject("GeminioApp", DefectRecombination);

template<>
InputParameters validParams<DefectRecombination>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("OtherDefect", "The other type of defect recombining with this defect");
  params.addRequiredCoupledVar("Recombination", "The recombination rate of both defects together");
  return params;
}


DefectRecombination::DefectRecombination(const
                                   InputParameters & parameters)
  :Kernel(parameters),
   _other_defect(coupledValue("OtherDefect")),
   _other_defect_var(coupled("OtherDefect")),
   _recombination_rate(coupledValue("Recombination"))
{}

Real
DefectRecombination::computeQpResidual()
{
  // Positive sign because negative source from weak form PDE
  return _test[_i][_qp] * _recombination_rate[_qp] * _other_defect[_qp] * _u[_qp];
}

Real
DefectRecombination::computeQpJacobian()
{
  // Positive sign because negative source from weak form PDE
  return _test[_i][_qp] * _recombination_rate[_qp] * _other_defect[_qp] * _phi[_j][_qp];
}

Real
DefectRecombination::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Positive sign because negative source from weak form PDE
  if (jvar == _other_defect_var)
    return _test[_i][_qp] * _recombination_rate[_qp] * _phi[_j][_qp] * _u[_qp];

  return 0.0;
}

#include "CoeffDiffusion.h"


registerMooseObject("GeminioApp", CoeffDiffusion);

template<>
InputParameters validParams<CoeffDiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addRequiredParam<Real>("diffusivity","diffusivity of this species");
  return params;
}

CoeffDiffusion::CoeffDiffusion(const InputParameters & parameters)
  : Diffusion(parameters),
    _diffusivity(getParam<Real>("diffusivity"))
{
}

Real
CoeffDiffusion::computeQpResidual()
{
  // Also... we're reusing the Diffusion Kernel's residual
  // so that we don't have to recode that.
  //  if (_u[_qp]>=0.0)
  return _diffusivity * Diffusion::computeQpResidual();
}

Real
CoeffDiffusion::computeQpJacobian()
{
  // Also... we're reusing the Diffusion Kernel's residual
  // so that we don't have to recode that.
  return _diffusivity * Diffusion::computeQpJacobian();
}

Real CoeffDiffusion::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
  return 0.0;
}

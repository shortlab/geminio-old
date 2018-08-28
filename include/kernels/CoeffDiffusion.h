#include "Diffusion.h"

#ifndef COEFFDIFFUSION_H
#define COEFFDIFFUSION_H

//Forward Declarations
class CoeffDiffusion;

template<>
InputParameters validParams<CoeffDiffusion>();

class CoeffDiffusion : public Diffusion
{
public:
  CoeffDiffusion(const  InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Material property of dispersion-diffusion coefficient.
  Real _diffusivity;

  // const MaterialProperty<Real> & _diffusivity;
};

#endif // COEFFDIFFUSION_H

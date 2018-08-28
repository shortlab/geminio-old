/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/

#ifndef GMATERIALCONSTANTS_H
#define GMATERIALCONSTANTS_H

#include "MaterialConstants.h"

class GMaterialConstants : public MaterialConstants
{
public:
  GMaterialConstants(const InputParameters & parameters);

  virtual Real absorbVV(int, int, int, Real) const { return 0.0; };
  virtual Real absorbVI(int, int, int, Real) const { return 0.0; };
  virtual Real absorbII(int, int, int, Real) const { return 0.0; };

  Real getAtomicVol() const { return _atomic_vol; };
  
protected:
  const Real _rho_d;
  const Real _i_bias;
  const Real _v_bias;

  const Real _atomic_vol;
};

template<>
InputParameters validParams<GMaterialConstants>();

#endif // GMATERIALCONSTANTS_H

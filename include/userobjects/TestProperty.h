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

#ifndef TESTPROPERTY_H
#define TESTPROPERTY_H

#include "GMaterialConstants.h"
#include "GeneralUserObject.h"

class TestProperty : public GMaterialConstants {
public:
  TestProperty(const InputParameters &parameters);

  Real absorb(int, int, MaterialParameters::Species, MaterialParameters::Species, Real, int, int) const;
  Real emit(int, int, Real, MaterialParameters::Species, MaterialParameters::Species, int, int) const;
  Real disl_ksq(int, MaterialParameters::Species, Real, bool = true) const;
  Real energy(int, MaterialParameters::Species, MaterialParameters::EType) const;
  Real D_prefactor(int, MaterialParameters::Species) const;
  Real diff(int, MaterialParameters::Species, Real) const;
  Real Ebinding(Real, MaterialParameters::Species, Real) const;

protected:
  const Real _rho_d;
  const Real _i_bias;
  const Real _v_bias;
  
  const Real _scale;
};

template <> InputParameters validParams<TestProperty>();

#endif

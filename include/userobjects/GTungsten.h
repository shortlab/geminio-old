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

#ifndef GTUNGSTEN_H
#define GTUNGSTEN_H

#include "GeneralUserObject.h"
#include "GMaterialConstants.h"

class GTungsten : public GMaterialConstants
{
public:
  GTungsten(const InputParameters & parameters);

  Real absorbVV(int, int, int, Real) const;
  Real absorbVI(int, int, int, Real) const;
  Real absorbII(int, int, int, Real) const;

  Real absorb(int, int, MaterialParameters::Species,MaterialParameters::Species, Real, int, int) const;
  Real emit(int, int, Real, MaterialParameters::Species, MaterialParameters::Species, int, int) const;
  Real diff(int, MaterialParameters::Species, Real) const;
  Real disl_ksq(int, MaterialParameters::Species, Real, bool mobile = true) const;

protected:  
  Real energy(int, MaterialParameters::Species, MaterialParameters::EType) const;
  Real D_prefactor(int, MaterialParameters::Species) const;

  const Real _burgers;
  const Real _rvi;
  const Real _Ev_formation;
  const Real _Ei_formation;
  const Real _Evb2;
  const Real _Eib2;
  const Real _Ei_binding_factor;
  const Real _Ev_binding_factor;
};

template<>
InputParameters validParams<GTungsten>();

#endif // GTUNGSTEN_H

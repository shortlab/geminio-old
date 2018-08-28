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

#ifndef GIRON_H
#define GIRON_H

#include "GeneralUserObject.h"
#include "GMaterialConstants.h"

/**
 * pure iron under neutron irradiation
 * calculate: Fe irradiation,  Neutron-induced swelling and embrittlement of pure
 * iron and pure nickel irradiated in the BN-350 and BOR-60 fast reactors
 * parameters: Efficient simulation of kinetics of radiation induced defects: A cluster dynamics approach
 */
class GIron : public GMaterialConstants
{
public:
  GIron(const InputParameters & parameters);

  Real absorb(int, int, MaterialParameters::Species, MaterialParameters::Species, Real, int, int) const;
  Real emit(int, int, Real, MaterialParameters::Species, MaterialParameters::Species, int, int) const;
  Real disl_ksq(int, MaterialParameters::Species, Real, bool mobile = true) const;
  Real energy(int, MaterialParameters::Species, MaterialParameters::EType) const;
  Real D_prefactor(int, MaterialParameters::Species) const;
  Real diff(int, MaterialParameters::Species, Real) const;
};

template<>
InputParameters validParams<GIron>();

#endif // GIRON_H

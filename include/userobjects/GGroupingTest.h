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

#ifndef GGROUPINGTEST_H
#define GGROUPINGTEST_H

#include "GeneralUserObject.h"
#include "GMaterialConstants.h"

/**
 * Based on test case (table 2) in GMIC++: Grouping method in C++: an efficient
 * method to solve large number of Master equations
 */
class GGroupingTest : public GMaterialConstants
{
public:
  GGroupingTest(const InputParameters & parameters);

  Real absorbVV(int, int, int, Real) const override;
  Real absorbVI(int, int, int, Real) const override;
  Real absorbII(int, int, int, Real) const override;

  Real absorb(int, int, MaterialParameters::Species, MaterialParameters::Species, Real, int, int) const override;
  Real emit(int, int, Real, MaterialParameters::Species, MaterialParameters::Species, int, int) const override;
  Real diff(int,  MaterialParameters::Species, Real) const override;
  Real disl_ksq(int, MaterialParameters::Species, Real, bool mobile = true) const override;
  
protected:
  Real energy(int, MaterialParameters::Species, MaterialParameters::EType) const;
  Real D_prefactor(int) const;
};

template<>
InputParameters validParams<GGroupingTest>();

#endif // GGROUPINGTEST_H

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

#ifndef GROUPINGTEST_H
#define GROUPINGTEST_H

#include "GeneralUserObject.h"
#include "GMaterialConstants.h"

/**
 * Based on test case (table 2) in GMIC++: Grouping method in C++: an efficient
 * method to solve large number of Master equations
 */
class GroupingTest : public GMaterialConstants
{
public:
  GroupingTest(const InputParameters & parameters);

  Real absorb(int, int, MaterialParameters::Species, MaterialParameters::Species, Real, int, int) const override;
  Real emit(int, int, Real, MaterialParameters::Species, MaterialParameters::Species, int, int) const override;
  Real disl_ksq(int, MaterialParameters::Species, Real, bool) const override { return 0.0; };
  Real diff(int,  MaterialParameters::Species,Real) const override;

protected:
  Real energy(int, MaterialParameters::Species, MaterialParameters::EType) const;
  Real D_prefactor(int, MaterialParameters::Species) const;

  Real _i_bias;
  Real _v_bias;
};

template<>
InputParameters validParams<GroupingTest>();

#endif // GROUPINGTEST_H

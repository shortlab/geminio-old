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

#ifndef MATERIALCONSTANTS_H
#define MATERIALCONSTANTS_H

#include "GeneralUserObject.h"
#include "MaterialParameters.h"

class MaterialConstants : public GeneralUserObject
{
public:
  MaterialConstants(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute() {}
  virtual void finalize() {}

  virtual Real absorb(int, int, MaterialParameters::Species, MaterialParameters::Species, Real, int, int) const = 0;
  virtual Real emit(int, int, Real, MaterialParameters::Species, MaterialParameters::Species, int, int) const = 0;
  virtual Real disl_ksq(int, MaterialParameters::Species, Real, bool) const = 0;
  virtual Real diff(int, MaterialParameters::Species, Real) const = 0;

protected:
  const Real _kB;
};

template<>
InputParameters validParams<MaterialConstants>();

#endif // MATERIALCONSTANTS_H

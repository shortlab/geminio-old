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

  virtual Real absorbVV(int,int,int,Real) const;
  virtual Real absorbVI(int,int,int,Real) const;
  virtual Real absorbII(int,int,int,Real) const;
  Real atomic_vol;

protected:
  Real _rho_d;
  Real _i_bias;
  Real _v_bias;
};

template<>
InputParameters validParams<GMaterialConstants>();

#endif // GMATERIALCONSTANTS_H

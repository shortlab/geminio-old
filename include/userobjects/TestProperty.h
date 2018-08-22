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

#include "GeneralUserObject.h"
#include "GMaterialConstants.h"

class TestProperty : public GMaterialConstants
{
public:
  TestProperty(const InputParameters & parameters);

  ~TestProperty(){}
  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  Real absorb(int,int,std::string,std::string,Real,int,int) const;
  Real emit(int,int,Real,std::string,std::string,int,int) const;
  Real disl_ksq(int,std::string,Real,int=1) const;
  Real energy(int,std::string,std::string) const;
  Real D_prefactor(int,std::string) const;
  Real diff(int, std::string,Real) const;
  Real Ebinding(Real,const char*,Real=1) const;

private:
  Real _rho_d;
  Real _i_bias;
  Real _v_bias;
};

template<>
InputParameters validParams<TestProperty>();

#endif

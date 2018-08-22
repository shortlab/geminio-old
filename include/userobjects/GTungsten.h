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

#ifndef GTungsten_H
#define GTungsten_H

#include "GeneralUserObject.h"
#include "GMaterialConstants.h"

class GTungsten : public GMaterialConstants
{
public:
  GTungsten(const InputParameters & parameters);

  ~GTungsten(){}
  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  Real absorb(int,int,std::string,std::string,Real,int,int) const;
  Real absorbVV(int,int,int,Real) const;
  Real absorbVI(int,int,int,Real) const;
  Real absorbII(int,int,int,Real) const;
  Real emit(int,int,Real,std::string,std::string,int,int) const;
  Real disl_ksq(int,std::string,Real,int=1) const;
  Real energy(int,std::string,std::string) const;
  Real D_prefactor(int,std::string="") const;
  Real diff(int, std::string,Real) const;

private:
  Real Ev_formation,Ei_formation,Evb2,Eib2;//description in cpp file
  Real Ei_binding_factor;
  Real Ev_binding_factor;
};

template<>
InputParameters validParams<GTungsten>();

#endif

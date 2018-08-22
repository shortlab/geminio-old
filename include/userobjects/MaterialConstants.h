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

class MaterialConstants : public GeneralUserObject
{
public:
  MaterialConstants(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  virtual Real absorb(int,int,std::string,std::string,Real,int,int) const;
  virtual Real emit(int,int,Real,std::string,std::string,int,int) const;
  virtual Real disl_ksq(int,std::string,Real,int=1) const;//1 denotes mobile
  virtual Real diff(int,std::string,Real) const;
};

template<>
InputParameters validParams<MaterialConstants>();

#endif // MATERIALCONSTANTS_H

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
#include "MooseMesh.h"
#include "MaterialConstants.h"

template<>
InputParameters validParams<MaterialConstants>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addClassDescription( "Object prototype to calculate material constants");
  return params;
}

MaterialConstants::MaterialConstants(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _kB(MaterialParameters::kB)
{
}

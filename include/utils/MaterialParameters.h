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

#ifndef MATERIALPARAMETERS_H
#define MATERIALPARAMETERS_H

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstring>

#include "Moose.h"

namespace MaterialParameters
{

enum class Species { V, I };
enum class EType { MIGRATION, BINDING };
  
Real absorb(int, int, Species, Species, Real, int, int);
Real emit(int, int, Real, Species, Species, int, int);
Real energy(int, Species, EType);
Real D_prefactor(int, Species);
Real diff(int, Species, Real);

}

#endif // MATERIALPARAMETERS_H

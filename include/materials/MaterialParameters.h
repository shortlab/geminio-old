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

Real absorb(int,int,std::string,std::string,Real,int,int);
Real emit(int,int,Real,std::string,std::string,int,int);
Real energy(int,std::string,std::string);
Real D_prefactor(int,std::string);
Real diff(int, std::string,Real);

#endif // MATERIALPARAMETERS_H

/**************************************
//Based on Table 1 in T. Jourdan's paper: Efficient simulation of kinetics of radiation induced defects: A cluster dynamics approach
*************************************/

#include "MaterialParameters.h"
#include "MooseError.h"
#include "libmesh/libmesh.h"

namespace MaterialParameters
{

const Real INF = 100.0;
const Real SCALE = 1.0;

//iron atom volume um^3
const Real Vatom = 1.182e-11;

Real
energy(int s, Species species, EType e_type)
{
  //unit:eV
  Real E = 0.0;
  if (species == Species::V && e_type == EType::MIGRATION)
  {
    ///* from literature
    switch (s)
    {
      case 1:
        E = 0.83;
        break;

      case 2:
        E = 0.62;
        break;

      case 3:
        E = 0.35;
        break;

      case 4:
        E = 0.48;
        break;

      default:
        E = INF;
    }
  }
  else if (species == Species::I && e_type == EType::MIGRATION)
  {
    switch (s)
    {
      case 1:
        E = 0.34;
        break;

      case 2:
        E = 0.42;
        break;

      case 3:
        E = 0.43;
        break;

      default:
        E = INF;
    }
  }
  else if (species == Species::V && e_type == EType::BINDING)
  {
    switch (s)
    {
      case 1:
        E = INF;
        break;

      case 2:
        E = 0.30;
        break;

      case 3:
        E = 0.37;
        break;

      case 4:
        E = 0.62;
        break;

      case 5:
        E = 0.73;
        break;

      default:
        // capillary law
        E = 2.2 + (0.3 - 2.2) / (std::pow(2.0, 2.0/3.0) - 1.0) * (std::pow(s, 2.0/3) - std::pow(s - 1.0, 2.0/3.0));
    }
  }
  else if (species == Species::I && e_type == EType::BINDING)
  {
    switch (s)
    {
      case 1:
        E = INF;
        break;

      case 2:
        E = 0.83;
        break;

      case 3:
        E = 0.92;
        break;

      case 4:
        E = 1.64;
        break;

      default:
        //E = 3.8 - 5.06*(pow(s,2.0/3)-pow(s-1,2.0/3));//capillary law
        // capillary law
        E = 3.64 - 4.78378 * (std::pow(s, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
    }
  }
  else
    mooseError("Energy not defined");

  return E;
}

Real
D_prefactor(int /*s*/, Species /*species*/)
{
  Real D0 = 8.2e5; // um^2/s
  return D0;
}

//size S1 and S2
Real
absorb(int S1, int S2, Species C1, Species C2, Real T, int tag1, int tag2)
{
  //tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
  if (tag1 == 0 && tag2 == 0)
    return 0.0;

  // recombination radius in um
  Real r_vi = 0.65e-3;

  Real r1 = std::pow(S1 * Vatom * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); //cluster effective radius
  Real r2 = std::pow(S2 * Vatom * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); //cluster effective radius
  Real D_s1 = D_prefactor(S1,C1) * std::exp(-energy(S1, C1, EType::MIGRATION) / (kB * T));
  Real D_s2 = D_prefactor(S2,C2) * std::exp(-energy(S2, C2, EType::MIGRATION) / (kB * T));
  return 4.0 * libMesh::pi * (D_s1 * tag1 + D_s2 * tag2) * (r1 + r2 + r_vi);
}

Real
diff(int S1, Species C1, Real T)
{
  //in um^2/s
  return D_prefactor(S1,C1) * std::exp(-energy(S1, C1, EType::MIGRATION) / (kB * T));
}

Real
emit(int S1, int S2, Real T, Species C1, Species /*C2*/, int tag1, int tag2)
{
  //for now only consider self species emission, S1 emits S2, S1==1
  if (S1 > S2 && S2 == 1)
    // unit: 1/s
    return  absorb(S1, S2, C1, C1, T, tag1, tag2) / (Vatom * std::pow(SCALE, 3)) * std::exp(-energy(S1, C1, EType::BINDING) / (kB * T));

  return 0.0;
}

}

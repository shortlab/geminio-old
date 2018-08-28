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
#include "BCCIronProperty.h"

registerMooseObject("GeminioApp", BCCIronProperty);

template<>
InputParameters validParams<BCCIronProperty>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addRequiredParam<Real>("dislocation", "dislocation density in #/um^2");
  params.addParam<Real>("i_disl_bias", 1.1, "dislocation bias for intersitials");
  params.addParam<Real>("v_disl_bias", 1.0, "dislocation bias for vacancies");

  // bcc iron atom volume um^3
  params.set<Real>("atomic_vol") = 1.181e-11;

  params.addClassDescription("Calculate BCC iron specific material properties");
  return params;
}

BCCIronProperty::BCCIronProperty(const InputParameters & parameters)
  : GMaterialConstants(parameters),
    _rho_d(getParam<Real>("dislocation")),
    _i_bias(getParam<Real>("i_disl_bias")),
    _v_bias(getParam<Real>("v_disl_bias")),
    _scale(1.0)
{
  _console << "BCCIronProperty constructed\n";
}

Real BCCIronProperty::energy(int s, MaterialParameters::Species species, MaterialParameters::EType e_type) const
{
  const Real INF = 100.0;

  // unit:eV
  if (e_type == MaterialParameters::EType::MIGRATION)
  {
    if (species == MaterialParameters::Species::V)
    {
      // from literature
      switch (s)
      {
        case 1:
          return 0.83;

        case 2:
          return 0.62;

        case 3:
          return 0.35;

        case 4:
          return 0.48;

        default:
          return INF;
      }
      // for test
      // if(s == 1) return 0.9; else return INF;
    }

    if (species == MaterialParameters::Species::I)
    {
      // data from literature
      switch (s)
      {
        case 1:
          return 0.34;

        case 2:
          return 0.42;

        case 3:
          return 0.43;

        default:
          return INF;
      }
      // for test
      // if( s<= 20) return 0.1+(s-1)/19.0*0.7; else return 1.1;//0.1+(s-1)/19.0*0.7;
    }
  }

  if (e_type == MaterialParameters::EType::BINDING)
  {
    if (species == MaterialParameters::Species::V)
    {
      switch (s){
        case 1:
          mooseWarning("Size should be larger than 1 to have a binding energy.");
          return 0.0;

        case 2:
          return 0.30;

        case 3:
          return 0.37;

        case 4:
          return 0.62;

        default:
          // larger cluster sizes
          return 2.2 - 3.2346 * (std::pow(s, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
      }
    }

    if (species == MaterialParameters::Species::I)
    {
      switch (s)
      {
        case 1:
          mooseWarning("Size should be larger than 1 to have a binding energy.");
          return 0.0;

        case 2:
          return 0.80;

        case 3:
          return 0.92;

        case 4:
          return 1.64;

        default:
          // larger cluster sizes
          return 3.8 - 5.06 * (std::pow(s, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
      }
    }
  }

  return INF;
}

Real
BCCIronProperty::D_prefactor(int s, MaterialParameters::Species species) const
{
  Real D0 = 0.0;
  if (species == MaterialParameters::Species::V)
  {
    // from literature
    switch (s)
    {
      case 1:
        D0 = 7.9e-7;//m^2/s
        break;

      case 2:
        D0 = 3.5e-8;
        break;

      default:
        D0 = 0.0;
    }
    // for test
    // if (s == 1) D0 = 2.0e-7; else D0 = 0.0;
  }
  else if (species == MaterialParameters::Species::I)
  {
    // from literature
    switch (s)
    {
      case 1:
        D0 = 1.3e-8;
        break;

      case 2:
        D0 = 351.6e-8; // this number is a little weird
        break;

      case 3:
        D0 = 12.1e-8;
        break;

      case 4:
        D0 = 12.3e-8;
        break;

      default:
        D0 = 9.0e-7 * std::pow(s, -0.6);
    }
    // for test
    // if (s <= 20) D0 = 2.0e-7/s; else D0 = 2.0e-7*pow(s,-0.7);
  }
  else
    D0 = 0.0;

  // change m^2/s to um^2/s
  return D0 * 1.0e12 * _scale * _scale;
}

//size S1 and S2
Real
BCCIronProperty::absorb(int S1, int S2, MaterialParameters::Species C1, MaterialParameters::Species C2, Real T, int tag1, int tag2) const
{
  // tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
  if (tag1 == 0 && tag2 == 0)
    return 0.0;

  Real r_vi = 0.65e-3 * _scale;
  Real r1 = std::pow(S1 * _atomic_vol * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); // cluster effective radius
  Real r2 = std::pow(S2 * _atomic_vol * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); // cluster effective radius
  Real D_s1 = D_prefactor(S1, C1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION)/ (_kB * T));
  Real D_s2 = D_prefactor(S2, C2) * std::exp(-energy(S2, C2, MaterialParameters::EType::MIGRATION)/ (_kB * T));
  return 4.0 * libMesh::pi * (D_s1 * tag1 + D_s2 * tag2) * (r1 + r2 + r_vi);
}

Real
BCCIronProperty::diff(int S1, MaterialParameters::Species C1, Real T) const
{
  // in um^2/s
  return D_prefactor(S1, C1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION)/ (_kB * T));
}

Real
BCCIronProperty::emit(int S1, int S2, Real T, MaterialParameters::Species C1, MaterialParameters::Species /*C2*/, int tag1, int tag2) const
{
  // for now only consider self species emission, S1 emits S2, S1==1
  Real emit_c = 0.0;
  if (S1 > S2 && S2==1)
    emit_c = absorb(S1, S2, C1, C1, T, tag1, tag2) / (_atomic_vol * std::pow(_scale, 3.0)) * std::exp(-energy(S1, C1, MaterialParameters::EType::BINDING)/ (_kB * T)); //unit:/s
  return emit_c;
}

Real
BCCIronProperty::disl_ksq(int S1, MaterialParameters::Species C1, Real T, bool mobile) const
{
  if (!mobile)
    return 0.0;

  const Real bias = C1 == MaterialParameters::Species::V ? _v_bias : _i_bias;
  return diff(S1, C1, T) * _rho_d * bias;
}

Real
BCCIronProperty::Ebinding(Real large, MaterialParameters::Species species, Real /*small*/) const
{
  // binding energy of small cluster, size 1, i.e. point defects
  int s = large;
  if (species == MaterialParameters::Species::V)
  {
    // vacancy type
    switch (s)
    {
      case 1:
        return 99999; // a huge value => impossible for emission

      case 2:
        return 0.30;

      case 3:
        return 0.37;

      case 4:
        return 0.62;

      default:
        return 2.2 - 3.2346 * (std::pow(s, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
    }
  }
  else if (species == MaterialParameters::Species::I)
  {
    // interstitial type
    switch (s)
    {
      case 1:
        return 99999;

      case 2:
        return 0.80;

      case 3:
        return 0.92;

      case 4:
        return 1.64;

      default:
        return 3.8 - 5.06 * (std::pow(s, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
    }
  }

  mooseError("Wrong argument to retrieve binding energy");
}

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
#include "TestProperty.h"

registerMooseObject("GeminioApp", TestProperty);

template<>
InputParameters validParams<TestProperty>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addRequiredParam<Real>("dislocation", "dislocation density in #/um^2");
  params.addParam<Real>("i_disl_bias", 1.1, "dislocation bias for interstitials");
  params.addParam<Real>("v_disl_bias", 1.0, "dislocation bias for vacancies");
  params.addClassDescription( "Calculate specific material properties");
  
  // iron atom volume um^3
  params.set<Real>("atomic_vol") = 1.182e-11;

  return params;
}

TestProperty::TestProperty(const InputParameters & parameters)
  : GMaterialConstants(parameters),
    _rho_d(getParam<Real>("dislocation")),
    _i_bias(getParam<Real>("i_disl_bias")),
    _v_bias(getParam<Real>("v_disl_bias")),
    _scale(1.0)
{
  _console << "TestProperty constructed\n";
}

Real TestProperty::energy(int s, MaterialParameters::Species species, MaterialParameters::EType e_type) const
{
  //unit:eV
  const Real INF = 100.0;
  
  if ((species == MaterialParameters::Species::V) && (e_type == MaterialParameters::EType::MIGRATION))
  {
    ///* from literature
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
  }
  else if ((species == MaterialParameters::Species::I) && (e_type == MaterialParameters::EType::MIGRATION))
  {
    switch (s)
    {
      case 1:
        return 0.34;
        break;

      case 2:
        return 0.42;
        break;

      case 3:
        return 0.43;
        break;

      default:
        return INF;
    }
  }
  else if ((species == MaterialParameters::Species::V) && (e_type == MaterialParameters::EType::BINDING))
  {
    switch (s){
      case 1:
        return INF;

      case 2:
        return 0.30;

      case 3:
        return 0.37;

      case 4:
        return 0.62;

      case 5:
        return 0.73;

      default:
        return 2.2 + (0.3-2.2)/(std::pow(2.0, 2.0/3.0)-1) * (std::pow(s, 2.0/3.0) - std::pow(s-1, 2.0/3.0));//capillary law
    }
  }
  else if ((species == MaterialParameters::Species::I) && (e_type == MaterialParameters::EType::BINDING))
  {
    switch (s)
    {
      case 1:
        return INF;

      case 2:
        return 0.83;

      case 3:
        return 0.92;

      case 4:
        return 1.64;

      default:
        // return 3.8 - 5.06*(std::pow(s, 2.0/3.0) - std::pow(s-1, 2.0/3.0));//capillary law
        // capillary law
        return 3.64 - 4.78378 * (std::pow(s, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
    }
  }
  else
    mooseError("Energy not defined");
}

Real
TestProperty::D_prefactor(int /*s*/, MaterialParameters::Species /*species*/) const
{
  Real D0 = 8.2e5; //um^2/s
  return D0;
}

//size S1 and S2
Real
TestProperty::absorb(int S1, int S2, MaterialParameters::Species C1, MaterialParameters::Species C2,Real T, int tag1, int tag2) const
{
  // tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
  if (tag1 == 0 && tag2 == 0)
    return 0.0;

  // recombination radius in um
  Real r_vi = 0.65e-3;
  Real r1 = std::pow(S1 * _atomic_vol * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); // cluster effective radius
  Real r2 = std::pow(S2 * _atomic_vol * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); // cluster effective radius
  Real D_s1 = D_prefactor(S1,C1) * std::exp(-energy(S1,C1,MaterialParameters::EType::MIGRATION) / (_kB * T));
  Real D_s2 = D_prefactor(S2,C2) * std::exp(-energy(S2,C2,MaterialParameters::EType::MIGRATION) / (_kB * T));
  return 4.0 * libMesh::pi * (D_s1 * tag1 + D_s2 * tag2) * (r1 + r2 + r_vi);
}

Real
TestProperty::diff(int S1, MaterialParameters::Species C1,Real T) const
{
  // in um^2/s
	return D_prefactor(S1,C1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION) / (_kB * T));
}

Real
TestProperty::emit(int S1, int S2, Real T, MaterialParameters::Species C1, MaterialParameters::Species /*C2*/, int tag1, int tag2) const
{
  // for now only consider self species emission, S1 emits S2, S1==1
  Real emit_c = 0.0;
  if (S1 > S2 && S2==1)
    emit_c = absorb(S1, S2, C1, C1, T, tag1, tag2)/(_atomic_vol * std::pow(_scale, 3.0)) * std::exp(-energy(S1, C1, MaterialParameters::EType::BINDING) / (_kB * T)); // unit:/s
  return emit_c;
}

Real
TestProperty::disl_ksq(int S1, MaterialParameters::Species C1, Real T, bool mobile) const
{
  if (!mobile)
    return 0.0;
  
  // dislocation sink rate k^2*Cj*Dj, return k^2*Dj in this function, where k^2= z*rho_d, P230/839 Was book
  const Real bias = C1 == MaterialParameters::Species::V ? _v_bias : _i_bias;
  return diff(S1, C1, T) * _rho_d * bias;
}

Real
TestProperty::Ebinding(Real large, MaterialParameters::Species species, Real /*small*/) const
{
  // binding energy of small cluster, size 1, i.e. point defects
  int s = large;
  if (species == MaterialParameters::Species::V)
  {
    // vacancy type
    switch (s)
    {
      case 1:
        // a huge value => impossible for emission
        return 99999; 

      case 2:
        return 0.30;

      case 3:
        return 0.37;

      case 4:
        return 0.62;

      default:
        return 2.2 - 3.2346 * (std::pow(s, 2.0/3.0) - std::pow(s - 1, 2.0/3.0));
    }
  }
  else if (species == MaterialParameters::Species::I)
  {
    switch (s)
    {
      // interstitial type
      case 1:
        return 99999;

      case 2:
        return 0.80;

      case 3:
        return 0.92;

      case 4:
        return 1.64;

      default:
        // return 3.8 - 5.06*(std::pow(s, 2.0/3.0) - std::pow(s-1, 2.0/3.0));
        // capillary law
        return 3.64 - 4.78378 * (std::pow(s, 2.0/3.0) - std::pow(s - 1, 2.0/3.0));
    }
  }
  else
    mooseError("Wrong argument to retrieve binding energy");
}

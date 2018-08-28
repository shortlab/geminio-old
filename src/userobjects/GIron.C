#include "MooseMesh.h"
#include "GIron.h"

registerMooseObject("GeminioApp", GIron);

template<>
InputParameters validParams<GIron>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addClassDescription( "Calculate material properties for pure iron under neutron irradiation");

  // iron atom volume um^3
  params.set<Real>("atomic_vol") = 1.205e-11;

  return params;
}

GIron::GIron(const InputParameters & parameters)
  : GMaterialConstants(parameters)
{
}

Real
GIron::energy(int s, MaterialParameters::Species species, MaterialParameters::EType e_type) const
{
  //unit:eV
  if (species == MaterialParameters::Species::V && e_type == MaterialParameters::EType::MIGRATION)
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
        return 100.0; // "INF"
    }
  }
  else if (species == MaterialParameters::Species::I && e_type == MaterialParameters::EType::MIGRATION)
  {
    switch (s)
    {
      case 1:
        return 0.34;

      case 2:
        return 0.42;

      case 3:
        return 0.43;

      default:
        return 100.0; // "INF"
    }
  }
  else if (species == MaterialParameters::Species::V && e_type == MaterialParameters::EType::BINDING)
  {
    switch (s)
    {
      case 1:
        return 100.0; // "INF"

      case 2:
        return 0.30;

      case 3:
        return 0.37;

      case 4:
        return 0.62;

      case 5:
        return 0.73;

      default:
        //capillary law
        return 2.2 + (0.3 - 2.2) / (std::pow(2.0, 2.0/3.0) - 1.0) * (std::pow(s, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
    }
  }
  else if (species == MaterialParameters::Species::I && e_type == MaterialParameters::EType::BINDING)
  {
    switch (s)
    {
      case 1:
        return 100.0; // "INF"

      case 2:
        return 0.83;

      case 3:
        return 0.92;

      case 4:
        return 1.64;

      default:
        //E = 3.8 - 5.06*(pow(s,2.0/3)-pow(s-1,2.0/3));//capillary law
        // capillary law
        return 3.64 - 4.78378 * (std::pow(s, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
    }
  }

  mooseError("Energy not defined");
}

Real
GIron::D_prefactor(int /*s*/, MaterialParameters::Species /*species*/) const
{
  Real D0 = 8.2e5; // um^2/s
  return D0;
}

//size S1 and S2
Real
GIron::absorb(int S1, int S2, MaterialParameters::Species C1, MaterialParameters::Species C2, Real T, int tag1, int tag2) const
{
  // tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
  if (tag1 == 0 && tag2 == 0)
    return 0.0;

  Real r_vi = 0.65e-3;//recombination radius in um
  Real r1 = std::pow(S1 * _atomic_vol * 3.0/(4.0 * libMesh::pi), 1.0/3.0); //cluster effective radius
  Real r2 = std::pow(S2 * _atomic_vol * 3.0/(4.0 * libMesh::pi), 1.0/3.0); //cluster effective radius
  Real D_s1 = D_prefactor(S1, C1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION)/ (_kB * T));
  Real D_s2 = D_prefactor(S2, C2) * std::exp(-energy(S2, C2, MaterialParameters::EType::MIGRATION)/ (_kB * T));
  return 4*libMesh::pi*(D_s1*tag1+D_s2*tag2)*(r1+r2+r_vi);
}

Real
GIron::diff(int S1, MaterialParameters::Species C1,Real T) const
{
  // in um^2/s
  return D_prefactor(S1,C1) * std::exp(-energy(S1,C1,MaterialParameters::EType::MIGRATION)/ (_kB * T));
}

Real
GIron::emit(int S1, int S2, Real T, MaterialParameters::Species C1, MaterialParameters::Species /*C2*/, int tag1, int tag2) const
{
  //for now only consider self species emission, S1 emits S2, S1==1
  Real emit_c = 0.0;
  if (S1 > S2 && S2 == 1)
    emit_c = absorb(S1, S2, C1, C1, T, tag1, tag2) / _atomic_vol * std::exp(-energy(S1, C1, MaterialParameters::EType::BINDING)/ (_kB * T));

  // unit:/s only emit point defect of the same species
  return emit_c;
}

// dislocation sink rate k^2*Cj*Dj, return k^2*Dj in this function, where k^2= z*rho_d, P230/839 Was book
Real
GIron::disl_ksq(int S1, MaterialParameters::Species C1, Real T, bool mobile) const
{
  if (!mobile)
    return 0.0;

  const Real bias = C1 == MaterialParameters::Species::V ? _v_bias : _i_bias;
  return diff(S1, C1, T) * _rho_d * bias;
}

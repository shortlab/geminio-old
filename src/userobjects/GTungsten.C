#include "MooseMesh.h"
#include "GTungsten.h"

registerMooseObject("GeminioApp", GTungsten);

template<>
InputParameters validParams<GTungsten>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addClassDescription("Calculate tungsten specific material properties");

  // tungsten atom volume um^3
  params.set<Real>("atomic_vol") = 1.5825e-11;

  return params;
}

GTungsten::GTungsten(const InputParameters & parameters)
  : GMaterialConstants(parameters),
    _burgers(2.7366e-4), // burgers vector (um) (sqrt(3)/2 * a0)
    _rvi(0.65e-3),       // vacancy - intersitial reaction distance (um))
    _Ev_formation(3.23), // vacancy formation energy eV
    _Ei_formation(9.96), // interstitial formation energy eV
    _Evb2(-0.1),         // binding energy for vacancy cluster size 2
    _Eib2(2.12),         // binding energy for interstitial cluster size 2
    _Ei_binding_factor((_Eib2 - _Ei_formation) / (std::pow(2.0, 2.0/3.0) - 1.0)),
    _Ev_binding_factor((_Evb2 - _Ev_formation) / (std::pow(2.0, 2.0/3.0) - 1.0))
{
}

Real
GTungsten::energy(int s, MaterialParameters::Species species, MaterialParameters::EType e_type) const
{
  //unit:eV
  //ref: Microstructural evolution of irradiated tungsten: Ab initio parameterisation of an OKMC model
  if (e_type == MaterialParameters::EType::MIGRATION)
  {
    if (species == MaterialParameters::Species::V)
      return 1.66; // eV

    if (species == MaterialParameters::Species::I)
      return 0.013; // .33; // eV
  }

  if (e_type == MaterialParameters::EType::BINDING)
  {
    if (species == MaterialParameters::Species::V)
    {
      switch (s)
      {
        case 1:
          return 100; // "INF"

        case 2:
          return _Evb2;

        case 3:
          return 0.04;

        case 4:
          return 0.64;

        case 5:
          return 0.72;

        case 6:
          return 0.89;

        case 7:
          return 0.72;

        default:
          return _Ev_formation + _Ev_binding_factor * (std::pow(s * 1.0, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));//capillary law
          // {
          //   Real Ef=1.77,gamma=1.0/1.6022e-7;//ev/um^2
          //   Real r = std::pow(s*_atomic_vol*3/4/libMesh::pi,1.0/3); //cluster effective radius
          //   return Ef-2*gamma*_atomic_vol/r;
          // }
      }
    }

    if (species == MaterialParameters::Species::I)
    {
      switch (s)
      {
        case 1:
          return 100; // "INF"

        case 2:
          return _Eib2;

        case 3:
          return 3.02;

        case 4:
          return 3.6;

        case 5:
          return 3.98;

        case 6:
          return 4.27;

        case 7:
          return 5.39;

        default:
          // capillary law
          return _Ei_formation + _Ei_binding_factor * (std::pow(s * 1.0, 2.0/3.0) - std::pow(s - 1.0, 2.0/3.0));
      }
    }
  }

  mooseError("Energy not defined");
}

Real
GTungsten::D_prefactor(int n, MaterialParameters::Species species) const
{
  //ref: Microstructural evolution of irradiated tungsten: Ab initio parameterization of an OKMC model
  // um^2/s
  if (species == MaterialParameters::Species::V)
    return 6.0096 * std::pow(10.0, 8.0 - 3.0 * n); // um^2/s
  else
    return 1.0016e5 * std::pow(1.0 * n, -0.5);
}

//size S1 and S2
Real
GTungsten::absorb(int S1, int S2, MaterialParameters::Species C1, MaterialParameters::Species C2,Real T, int tag1, int tag2) const
{
  // tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
  if (tag1 == 0 && tag2 == 0)
    return 0.0;

  int S = S1 > S2 ? S1 : S2;
  Real w = std::pow(48.0 * libMesh::pi * libMesh::pi /_atomic_vol/_atomic_vol * S, 1.0/3.0);
  Real D_s1 = D_prefactor(S1, C1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION) / (_kB * T));
  Real D_s2 = D_prefactor(S2, C2) * std::exp(-energy(S2, C2, MaterialParameters::EType::MIGRATION) / (_kB * T));
  return w*_atomic_vol*(D_s1*tag1+D_s2*tag2);//add _atomic_vol for unit concern,  P5/19 in ref

  //test the other expression, larger than the one above (3 times)
  /*
    Real r_vi = 0.65e-3;//recombination radius in um
    Real r1 = pow(S1*_atomic_vol*3/4/libMesh::pi,1.0/3); //cluster effective radius
    Real r2 = pow(S2*_atomic_vol*3/4/libMesh::pi,1.0/3); //cluster effective radius
    return 4*libMesh::pi*(D_s1*tag1+D_s2*tag2)*(r1+r2+r_vi);
  */
}

//vv reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
Real
GTungsten::absorbVV(int S1, int S2, int flag, Real T) const
{
  Real result = 0.0;
  switch (flag)
  {
    case 1:
    {
      const Real w = std::pow(48.0 * libMesh::pi * libMesh::pi / _atomic_vol / _atomic_vol * S2, 1.0/3);
      const Real D_s1 = D_prefactor(S1,MaterialParameters::Species::V) *
                        std::exp(-energy(S1, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s1;
      break;
    }
    case 2:
    {
      const Real w = std::pow(48.0 * libMesh::pi * libMesh::pi / _atomic_vol / _atomic_vol * S1, 1.0/3.0);
      const Real D_s2 = D_prefactor(S2, MaterialParameters::Species::V) *
                        std::exp(-energy(S2, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s2;
      break;
    }
    case 3:
    {
      int S = S1 > S2 ? S1 : S2;
      const Real w = std::pow(48.0 * libMesh::pi * libMesh::pi / _atomic_vol / _atomic_vol * S, 1.0/3.0);
      const Real D_s1 = D_prefactor(S1, MaterialParameters::Species::V) *
                        std::exp(-energy(S1, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      const Real D_s2 = D_prefactor(S2, MaterialParameters::Species::V) *
                        std::exp(-energy(S2, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      // add _atomic_vol for unit concern, P5/19 in ref
      result = w * _atomic_vol * (D_s1 + D_s2);
      break;
    }
  }

  return result;
}

//vi reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
Real
GTungsten::absorbVI(int S1, int S2, int flag, Real T) const
{
  Real result = 0.0;
  switch (flag)
  {
    case 1:
    {
      const Real w = std::pow(4*libMesh::pi/_atomic_vol/_burgers*S2,1.0/2);
      const Real D_s1 = D_prefactor(S1,MaterialParameters::Species::V) * std::exp(-energy(S1, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = _v_bias*w*_atomic_vol*D_s1;
      break;
    }
    case 2:
    {
      const Real w = std::pow(48.0 * libMesh::pi * libMesh::pi / _atomic_vol / _atomic_vol * S1, 1.0/3);
      const Real D_s2 = D_prefactor(S2, MaterialParameters::Species::I) *
                        std::exp(-energy(S2, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s2;
      break;
    }
    case 3:
    {
      const Real D_s1 = D_prefactor(S1, MaterialParameters::Species::V) *
                        std::exp(-energy(S1, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      const Real D_s2 = D_prefactor(S2, MaterialParameters::Species::I) *
                        std::exp(-energy(S2, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      //ref: Mean field rate theory and object kinetic monte carlo: a comparison of kinetic models
      result = 4.0 * libMesh::pi * _rvi / _atomic_vol * (D_s1 + D_s2);
      break;
    }
  }
  return result;
}

//ii reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
Real
GTungsten::absorbII(int S1, int S2, int flag, Real T) const
{
  Real result = 0.0;
  switch (flag)
  {
    case 1:
    {
      const Real w = std::pow(4.0 * libMesh::pi / _atomic_vol / _burgers * S2, 1.0/2.0);
      Real D_s1 = D_prefactor(S1,MaterialParameters::Species::I) *
                  std::exp(-energy(S1, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = _i_bias*w*_atomic_vol*D_s1;
      break;
    }
    case 2:
    {
      const Real w = std::pow(4.0 * libMesh::pi / _atomic_vol / _burgers * S1, 1.0/2.0);
      Real D_s2 = D_prefactor(S2,MaterialParameters::Species::I) *
                  std::exp(-energy(S2, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = _i_bias * w * _atomic_vol * D_s2;
      break;
    }
    case 3:
    {
      int S = (S1>S2)? S1: S2;
      const Real w = std::pow(4.0 * libMesh::pi / _atomic_vol / _burgers * S,1.0/2);
      Real D_s1 = D_prefactor(S1,MaterialParameters::Species::I) * std::exp(-energy(S1, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      Real D_s2 = D_prefactor(S2,MaterialParameters::Species::I) * std::exp(-energy(S2, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = _i_bias*w*_atomic_vol*(D_s1+D_s2);//ref: Mean field rate theory and object kinetic monte carlo: a comparison of kinetic models
      break;
    }
  }
  return result;
}

Real
GTungsten::diff(int S1, MaterialParameters::Species C1,Real T) const
{
  //in um^2/s
  return D_prefactor(S1,C1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION) / (_kB * T));
}

Real
GTungsten::emit(int S1, int S2, Real T, MaterialParameters::Species C1, MaterialParameters::Species /*C2*/, int tag1, int tag2) const
{
  //for now only consider self species emission, S1 emits S2, S1==1
  if (C1 == MaterialParameters::Species::I)
    //intersitial cluster doesn't emit.
    return 0.0;

  Real emit_c = 0.0;
  if (S1 > S2 && S2==1)
    emit_c = absorb(S1, S2, C1, C1, T, tag1, tag2) / _atomic_vol * std::exp(-energy(S1, C1, MaterialParameters::EType::BINDING) / (_kB * T));

  // unit: 1/s only emit point defect of the same species
  return emit_c;
}

Real
GTungsten::disl_ksq(int S1, MaterialParameters::Species C1, Real T, bool mobile) const
{
  if (!mobile)
    return 0.0;

  const Real bias = C1 == MaterialParameters::Species::V ? _v_bias : _i_bias;
  return diff(S1,C1,T) * _rho_d * bias;
}

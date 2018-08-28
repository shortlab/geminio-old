#include "MooseMesh.h"
#include "GGroupingTest.h"

registerMooseObject("GeminioApp", GGroupingTest);

template<>
InputParameters validParams<GGroupingTest>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addClassDescription("Calculate specific material properties");
  
  // iron atom volume um^3
  params.set<Real>("atomic_vol") = 1.205e-11;

  return params;
}

GGroupingTest::GGroupingTest(const InputParameters & parameters)
  : GMaterialConstants(parameters)
{
  // _console << "GGroupingTest constructed\n";
}

/*
Real GGroupingTest::energy(int s, MaterialParameters::Species species, MaterialParameters::EType e_type) const
{
  // unit:eV
  const Real INF = 100.0;
  if (species == MaterialParameters::Species::V && e_type == MaterialParameters::EType::MIGRATION)
  {
    // from literature
    switch(s)
    {
      case 1:
        // TODO: use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
        return 1.1;

      default:
        return INF;
    }
  }
  else if (species == MaterialParameters::Species::I && e_type == MaterialParameters::EType::MIGRATION)
  {
    return INF;
  }
  else if (species == MaterialParameters::Species::V && e_type == MaterialParameters::EType::BINDING)
  {
    Real Ef=1.77, gamma=1.0/1.6022e-7;//ev/um^2
    Real r = pow(s*_atomic_vol*3/4/libMesh::pi, 1.0/3); //cluster effective radius
    return Ef-2*gamma*_atomic_vol/r;
  }
  else if (species == MaterialParameters::Species::I && e_type == MaterialParameters::EType::BINDING)
  {
    return INF;
  }
  else
    mooseError("Energy not defined");
    
  return E;
}
*/

// test interstitial part
Real 
GGroupingTest::energy(int s, MaterialParameters::Species species, MaterialParameters::EType e_type) const
{
  //unit:eV
  const Real INF = 100.0;
  if (e_type == MaterialParameters::EType::MIGRATION)
  {
    if (species == MaterialParameters::Species::I)
    {
      // from literature
      switch (s)
      {
        case 1:
          // 0.34; TODO:use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
          return 1.1;

        case 2:
          return 0.42;

        case 3:
          return 0.43;

        default:
          return INF;
      }
    }
    else if (species == MaterialParameters::Species::V)
    {
      switch (s)
      {
        case 1:
          // TODO:use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
          return 1.1;

        case 2:
          return 0.62;

        case 3:
          return 0.35;

        default:
          return INF;
      }
    }
  }

  if (e_type == MaterialParameters::EType::BINDING)
  {
    if (species == MaterialParameters::Species::I) 
    {
      if (s == 1)
        return INF;

      Real Ef = 1.77, gamma = 1.0 / 1.6022e-7; // ev/um^2
      // cluster effective radius
      Real r = std::pow(s * _atomic_vol * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); 
      return Ef - 2 * gamma * _atomic_vol / r;
    }
    if (species == MaterialParameters::Species::V)
    {
      if (s == 1)
        return INF;
        
      Real Ef = 1.77, gamma = 1.0 / 1.6022e-7; // ev/um^2
      // cluster effective radius
      Real r = std::pow(s  *_atomic_vol * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); 
      return Ef-2*gamma*_atomic_vol/r;
    }
  }

  mooseError("Energy not defined");
}


Real
GGroupingTest::D_prefactor(int /*s*/) const
{
  Real D0 = 1.0e6; // um^2/s
  return D0;
}

//size S1 and S2
Real
GGroupingTest::absorb(int S1, int S2, MaterialParameters::Species C1, MaterialParameters::Species C2, Real T, int tag1, int tag2) const
{
  //tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
  if (tag1 == 0 && tag2 == 0)
    return 0.0;

  int S = (S1 > S2) ? S1 : S2;
  Real w = std::pow(48.0 * libMesh::pi * libMesh::pi / _atomic_vol / _atomic_vol * S, 1.0/3.0);
  Real D_s1 = D_prefactor(S1) * std::exp(-energy(S1,C1, MaterialParameters::EType::MIGRATION) / (_kB * T));
  Real D_s2 = D_prefactor(S2) * std::exp(-energy(S2,C2, MaterialParameters::EType::MIGRATION) / (_kB * T));

  // add _atomic_vol for unit concern,  P5/19 in ref
  return w * _atomic_vol * (D_s1 * tag1 + D_s2 * tag2);

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
GGroupingTest::absorbVV(int S1, int S2, int flag, Real T) const
{
  Real result = 0.0;
  int S = (S1 > S2) ? S1 : S2;
  Real w = std::pow(48.0 * libMesh::pi * libMesh::pi / _atomic_vol / _atomic_vol * S, 1.0/3.0);
  switch (flag)
  {
    case 1:
    {
      Real D_s1 = D_prefactor(S1) * std::exp(-energy(S1, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s1;
      break;
    }
    case 2:
    {
      Real D_s2 = D_prefactor(S2) * std::exp(-energy(S2, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s2;
      break;
    }
    case 3:
    {
      Real D_s1 = D_prefactor(S1) * std::exp(-energy(S1, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      Real D_s2 = D_prefactor(S2) * std::exp(-energy(S2, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w*_atomic_vol*(D_s1+D_s2);//add _atomic_vol for unit concern,  P5/19 in ref
      break;
    }
  }
  return result;
}

//vi reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
Real
GGroupingTest::absorbVI(int S1, int S2, int flag, Real T) const
{
  Real result = 0.0;
  int S = (S1>S2)? S1: S2;
  Real w = pow(48.0*libMesh::pi*libMesh::pi/_atomic_vol/_atomic_vol*S,1.0/3);
  switch (flag)
  {
    case 1:
    {
      Real D_s1 = D_prefactor(S1) * std::exp(-energy(S1, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s1;
      break;
    }
    case 2:
    {
      Real D_s2 = D_prefactor(S2) * std::exp(-energy(S2, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s2;
      break;
    }
    case 3:
    {
      Real D_s1 = D_prefactor(S1) * std::exp(-energy(S1, MaterialParameters::Species::V, MaterialParameters::EType::MIGRATION) / (_kB * T));
      Real D_s2 = D_prefactor(S2) * std::exp(-energy(S2, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w*_atomic_vol*(D_s1+D_s2);//ref: Mean field rate theory and object kinetic monte carlo: a comparison of kinetic models
      break;
    }
  }
  return result;
}

//ii reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
Real
GGroupingTest::absorbII(int S1, int S2, int flag, Real T) const
{
  Real result = 0.0;
  int S = (S1>S2)? S1: S2;
  Real w = pow(48.0*libMesh::pi*libMesh::pi/_atomic_vol/_atomic_vol*S,1.0/3);
  switch (flag)
  {
    case 1:
    {
      Real D_s1 = D_prefactor(S1) * std::exp(-energy(S1, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s1;
      break;
    }
    case 2:
    {
      Real D_s2 = D_prefactor(S2) * std::exp(-energy(S2, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      result = w * _atomic_vol * D_s2;
      break;
    }
    case 3:
    {
      Real D_s1 = D_prefactor(S1) * std::exp(-energy(S1, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      Real D_s2 = D_prefactor(S2) * std::exp(-energy(S2, MaterialParameters::Species::I, MaterialParameters::EType::MIGRATION) / (_kB * T));
      //ref: Mean field rate theory and object kinetic Monte Carlo: a comparison of kinetic models
      result = w * _atomic_vol * (D_s1 + D_s2);
      break;
    }
  }
  return result;
}

Real
GGroupingTest::diff(int S1, MaterialParameters::Species C1, Real T) const
{
  //in um^2/s
  return D_prefactor(S1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION) / (_kB * T));
}

Real
GGroupingTest::emit(int S1, int S2, Real T, MaterialParameters::Species C1, MaterialParameters::Species /*C2*/, int tag1, int tag2) const
{
  //for now only consider self species emission, S1 emits S2, S1==1
  Real emit_c = 0.0;
  if (S1 > S2 && S2==1)
    emit_c = absorb(S1, S2, C1, C1, T, tag1, tag2)/(_atomic_vol) * 
             std::exp(-energy(S1, C1, MaterialParameters::EType::BINDING) / (_kB * T));

  // unit:/s only emit point defect of the same species
  return emit_c;
}

Real 
GGroupingTest::disl_ksq(int S1, MaterialParameters::Species C1, Real T, bool mobile) const 
{
  if (!mobile)
    return 0.0;
    
  const Real bias = C1 == MaterialParameters::Species::V ? _v_bias : _i_bias;
  return diff(S1, C1, T) * _rho_d * bias;
}

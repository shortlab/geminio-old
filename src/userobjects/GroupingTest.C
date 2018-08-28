#include "MooseMesh.h"
#include "GroupingTest.h"

registerMooseObject("GeminioApp", GroupingTest);

template<>
InputParameters validParams<GroupingTest>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addParam<Real>("i_disl_bias", 1.1, "dislocation bias for interstitials");
  params.addParam<Real>("v_disl_bias", 1.0, "dislocation bias for vacancies");
    
  // iron atom volume um^3
  params.set<Real>("atomic_vol") = 1.205e-11;

  params.addClassDescription( "Calculate specific material properties");
  return params;
}

GroupingTest::GroupingTest(const InputParameters & parameters)
  : GMaterialConstants(parameters),
    _i_bias(getParam<Real>("i_disl_bias")),
    _v_bias(getParam<Real>("v_disl_bias"))
{
  _console << "GroupingTest constructed\n";
}

Real
GroupingTest::energy(int s, MaterialParameters::Species species, MaterialParameters::EType e_type) const
{
  //unit:eV
  const Real INF = 100.0;
  if (e_type == MaterialParameters::EType::MIGRATION)
  {
    if (species == MaterialParameters::Species::V)
    {
      // from literature
      switch (s)
      {
        case 1:
          // TODO: use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
          return 1.1;

        default:
          return INF;
      }
    }
    else if (species == MaterialParameters::Species::I)
    {
      switch (s)
      {
        case 1:
          // TODO: use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
          return 1.1;

        default:
          return INF;
      }
    }
  }
  
  if (e_type == MaterialParameters::EType::BINDING)
  {
    if (species == MaterialParameters::Species::V)
    {
      Real Ef=1.77, gamma = 1.0 / 1.6022e-7; // ev/um^2
      // cluster effective radius
      Real r = std::pow(s * _atomic_vol * 3.0 / (4.0 * libMesh::pi), 1.0/3.0); 
      return Ef - 2 * gamma * _atomic_vol / r;
    }
    else if (species == MaterialParameters::Species::I)
    {
      Real Ef = 1.77, gamma = 1.0 / 1.6022e-7; // ev/um^2
      // cluster effective radius
      Real r = std::pow(s * _atomic_vol * 3.0 / (4 * libMesh::pi), 1.0/3.0); 
      return Ef - 2 * gamma * _atomic_vol / r;
    }
  }

  mooseError("Energy not defined");
}

Real 
GroupingTest::D_prefactor(int /*s*/, MaterialParameters::Species /*species*/) const
{
  Real D0 = 1.0e6; // um^2/s
  return D0;
}

// size S1 and S2
Real 
GroupingTest::absorb(int S1, int S2, MaterialParameters::Species C1, MaterialParameters::Species C2,Real T, int tag1, int tag2) const
{
  // tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
  if (tag1 == 0 && tag2 == 0) 
    return 0.0;
    
  int S = S1 > S2 ? S1 : S2;
  Real w = std::pow(48.0 * libMesh::pi * libMesh::pi / _atomic_vol / _atomic_vol * S, 1.0/3.0);
  Real D_s1 = D_prefactor(S1, C1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION) / (_kB * T));
  Real D_s2 = D_prefactor(S2, C2) * std::exp(-energy(S2, C2, MaterialParameters::EType::MIGRATION) / (_kB * T));
  return w*_atomic_vol*(D_s1*tag1+D_s2*tag2);//add _atomic_vol for unit concern,  P5/19 in ref

//test the other expression, larger than the one above (3 times)
/*
    Real r_vi = 0.65e-3;//recombination radius in um
    Real r1 = pow(S1*_atomic_vol*3/4/libMesh::pi, 1.0/3); //cluster effective radius
    Real r2 = pow(S2*_atomic_vol*3/4/libMesh::pi, 1.0/3); //cluster effective radius
    return 4*libMesh::pi*(D_s1*tag1+D_s2*tag2)*(r1+r2+r_vi);
*/
}

Real 
GroupingTest::diff(int S1, MaterialParameters::Species C1, Real T) const 
{
  // in um^2/s
  return D_prefactor(S1, C1) * std::exp(-energy(S1, C1, MaterialParameters::EType::MIGRATION) / (_kB * T));
}

Real
GroupingTest::emit(int S1, int S2, Real T, MaterialParameters::Species C1, MaterialParameters::Species /*C2*/, int tag1, int tag2) const
{
  //for now only consider self species emission, S1 emits S2, S1==1
  Real emit_c = 0.0;
  if (S1 > S2 && S2 == 1)
    emit_c = absorb(S1, S2, C1, C1, T, tag1, tag2)/(_atomic_vol) * 
             std::exp(-energy(S1, C1, MaterialParameters::EType::BINDING) / (_kB * T));

  // unit:/s only emit point defect of the same species
  return emit_c;
}

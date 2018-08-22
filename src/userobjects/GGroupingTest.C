//Based on test case (table 2) in GMIC++: Grouping method in C++: an efficient method to solve large number of Master equations
#include "MooseMesh.h"
#include "GGroupingTest.h"

#define INF 100
#define SCALE 1 //change unit from um
#define PI 3.14159265359
#define Vatom 1.205e-11 //iron atom volume um^3
#define Boltz_const 8.6173315e-5 //boltzmann constant eV/K

template<>
InputParameters validParams<GGroupingTest>()
{
  InputParameters params = validParams<GMaterialConstants>();
  //params.addParam<Real>("i_disl_bias",1.1,"dislocation bias for interstitials");
  //params.addParam<Real>("v_disl_bias",1.0,"dislocation bias for vacancies");
  params.addClassDescription( "Calculate specific material properties");
  return params;
}

GGroupingTest::GGroupingTest(const InputParameters & parameters)
: GMaterialConstants(parameters)
{
//printf("GGroupingTest constructed\n");
  atomic_vol = Vatom;
}

void GGroupingTest::initialize()
{
}

void GGroupingTest::execute()
{
}

void GGroupingTest::finalize()
{
}
/*
Real GGroupingTest::energy(int s,std::string species, std::string Etype) const{//unit:eV
    Real E=0.0;
    if ((species == "V") && (Etype == "migration")){
//// from literature
        switch(s){
            case 1:
            {
                E = 1.1;//TODO:use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
                break;
            }
            default:
                E = INF;
        }
    }
    else if ((species == "I") && (Etype == "migration")){
        return INF;
    }
    else if ((species == "V") && (Etype == "binding")){
        Real Ef=1.77,gamma=1.0/1.6022e-7;//ev/um^2
        Real r = pow(s*Vatom*3/4/PI,1.0/3); //cluster effective radius
        return Ef-2*gamma*Vatom/r;
    }
    else if ((species == "I") && (Etype == "binding")) {
        return INF;
    }
    else
        mooseError("Energy not defined for " + Etype + " " + species);
    return E;
}
*/
//test intersititial part
Real GGroupingTest::energy(int s,std::string species, std::string Etype) const{//unit:eV
    Real E=0.0;
    if ((species == "I") && (Etype == "migration")){
//// from literature
        switch(s){
            case 1:
            {
                E = 1.1;//0.34;//TODO:use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
                break;
            }
            case 2:
            {
                E = 0.42;
                break;
            }
            case 3:
            {
                E = 0.43;
                break;
            }
            default:
                E = INF;
        }
    }
    else if ((species == "V") && (Etype == "migration")){
        switch(s){
            case 1:
            {
                E = 1.1;//1.1;//TODO:use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
                break;
            }
            case 2:
            {
                E = 0.62;
                break;
            }
            case 3:
            {
                E = 0.35;
                break;
            }
            default:
                E = INF;
        }
    }
    else if ((species == "I") && (Etype == "binding")){
      if (s == 1)
        return INF;
      Real Ef=1.77,gamma=1.0/1.6022e-7;//ev/um^2
      Real r = pow(s*Vatom*3/4/PI,1.0/3); //cluster effective radius
      return Ef-2*gamma*Vatom/r;
    }
    else if ((species == "V") && (Etype == "binding")) {
      if (s==1)
        return INF;
      Real Ef = 1.77, gamma = 1.0/1.6022e-7; // ev/um^2
      Real r = pow(s*Vatom*3/4/PI,1.0/3); // cluster effective radius
      return Ef-2*gamma*Vatom/r;
    }
    else
      mooseError("Energy not defined for " + Etype + " " + species);
    return E;
}


Real
GGroupingTest::D_prefactor(int /*s*/, std::string /*species*/) const
{
  Real D0 = 1.0e6;//um^2/s
  return D0;
}

//size S1 and S2
Real
GGroupingTest::absorb(int S1, int S2, std::string C1, std::string C2, Real T, int tag1, int tag2) const
{
  //tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
  if (tag1 == 0 && tag2 == 0)
    return 0.0;

  int S = (S1 > S2) ? S1 : S2;
  Real w = pow(48.0*PI*PI/Vatom/Vatom*S,1.0/3);
  Real D_s1 = D_prefactor(S1,C1) * exp(-energy(S1,C1, "migration")/Boltz_const/T);
  Real D_s2 = D_prefactor(S2,C2) * exp(-energy(S2,C2, "migration")/Boltz_const/T);

  // add Vatom for unit concern,  P5/19 in ref
  return w * Vatom * (D_s1 * tag1 + D_s2 * tag2);

  //test the other expression, larger than the one above (3 times)
  /*
      Real r_vi = 0.65e-3;//recombination radius in um
      Real r1 = pow(S1*Vatom*3/4/PI,1.0/3); //cluster effective radius
      Real r2 = pow(S2*Vatom*3/4/PI,1.0/3); //cluster effective radius
      return 4*PI*(D_s1*tag1+D_s2*tag2)*(r1+r2+r_vi);
  */
}

//vv reaction; flag=0: both immobile; flag=1: first mobile; flag=2: second mobile; flag=3: both mobile
Real
GGroupingTest::absorbVV(int S1, int S2, int flag, Real T) const
{
  Real result = 0.0;
  int S = (S1 > S2) ? S1 : S2;
  Real w = pow(48.0*PI*PI/Vatom/Vatom*S,1.0/3);
  switch (flag)
  {
    case 1:
    {
      Real D_s1 = D_prefactor(S1)*exp(-energy(S1,"V","migration")/Boltz_const/T);
      result = w*Vatom*D_s1;
      break;
    }
    case 2:
    {
      Real D_s2 = D_prefactor(S2)*exp(-energy(S2,"V","migration")/Boltz_const/T);
      result = w*Vatom*D_s2;
      break;
    }
    case 3:
    {
      Real D_s1 = D_prefactor(S1)*exp(-energy(S1,"V","migration")/Boltz_const/T);
      Real D_s2 = D_prefactor(S2)*exp(-energy(S2,"V","migration")/Boltz_const/T);
      result = w*Vatom*(D_s1+D_s2);//add Vatom for unit concern,  P5/19 in ref
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
  Real w = pow(48.0*PI*PI/Vatom/Vatom*S,1.0/3);
  switch (flag)
  {
    case 1:
    {
      Real D_s1 = D_prefactor(S1)*exp(-energy(S1,"V","migration")/Boltz_const/T);
      result = w*Vatom*D_s1;
      break;
    }
    case 2:
    {
      Real D_s2 = D_prefactor(S2)*exp(-energy(S2,"I","migration")/Boltz_const/T);
      result = w*Vatom*D_s2;
      break;
    }
    case 3:
    {
      Real D_s1 = D_prefactor(S1)*exp(-energy(S1,"V","migration")/Boltz_const/T);
      Real D_s2 = D_prefactor(S2)*exp(-energy(S2,"I","migration")/Boltz_const/T);
      result = w*Vatom*(D_s1+D_s2);//ref: Mean field rate theory and object kinetic monte carlo: a comparison of kinetic models
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
  Real w = pow(48.0*PI*PI/Vatom/Vatom*S,1.0/3);
  switch (flag)
  {
    case 1:
    {
      Real D_s1 = D_prefactor(S1)*exp(-energy(S1,"I","migration")/Boltz_const/T);
      result = w*Vatom*D_s1;
      break;
    }
    case 2:
    {
      Real D_s2 = D_prefactor(S2)*exp(-energy(S2,"I","migration")/Boltz_const/T);
      result = w*Vatom*D_s2;
      break;
    }
    case 3:
    {
      Real D_s1 = D_prefactor(S1)*exp(-energy(S1,"I","migration")/Boltz_const/T);
      Real D_s2 = D_prefactor(S2)*exp(-energy(S2,"I","migration")/Boltz_const/T);
      result = w*Vatom*(D_s1+D_s2);//ref: Mean field rate theory and object kinetic monte carlo: a comparison of kinetic models
      break;
    }
  }
  return result;
}

Real
GGroupingTest::diff(int S1, std::string C1,Real T) const
{
  //in um^2/s
  return D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
}

Real
GGroupingTest::emit(int S1, int S2, Real T, std::string C1, std::string /*C2*/, int tag1, int tag2) const
{
  //for now only consider self species emission, S1 emits S2, S1==1
  Real emit_c = 0.0;
  if (S1 > S2 && S2==1)
    emit_c = absorb(S1,S2,C1,C1,T,tag1,tag2)/(Vatom) *exp(-energy(S1,C1,"binding")/Boltz_const/T);

  // unit:/s only emit point defect of the same species
  return emit_c;
}

Real GGroupingTest::disl_ksq(int S1, std::string C1, Real T, int tag) const {
   Real bias = (! C1.compare("V"))? _v_bias : _i_bias;
   return tag * diff(S1,C1,T) * _rho_d * bias;
}

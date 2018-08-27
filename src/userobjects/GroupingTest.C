//Based on test case (table 2) in GMIC++: Grouping method in C++: an efficient method to solve large number of Master equations
#include "MooseMesh.h"
#include "GroupingTest.h"

#define INF 100
#define SCALE 1 //change unit from um
#define PI 3.14159265359
#define Vatom 1.205e-11 //iron atom volume um^3
#define Boltz_const 8.6173315e-5 //boltzmann constant eV/K

registerMooseObject("GeminioApp", GroupingTest);

template<>
InputParameters validParams<GroupingTest>()
{
  InputParameters params = validParams<GMaterialConstants>();
  params.addParam<Real>("i_disl_bias",1.1,"dislocation bias for intersitials");
  params.addParam<Real>("v_disl_bias",1.0,"dislocation bias for vacancies");
  params.addClassDescription( "Calculate specific material properties");
  return params;
}

GroupingTest::GroupingTest(const InputParameters & parameters)
: GMaterialConstants(parameters),
 _i_bias(getParam<Real>("i_disl_bias")),
 _v_bias(getParam<Real>("v_disl_bias"))
{
printf("GroupingTest constructed\n");
}

Real
GroupingTest::energy(int s,std::string species, std::string Etype) const
{
  //unit:eV
  Real E=0.0;
  if ((species == "V") && (Etype == "migration"))
  {
    // from literature
    switch (s)
    {
      case 1:
        E = 1.1;//TODO:use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
        break;

      default:
        E = INF;
    }
  }
  else if ((species == "I") && (Etype == "migration")){
    switch (s)
    {
      case 1:
        E = 1.1;//TODO:use 1.3 from another ref: Cluster dynamics simulation of point defect clusters in neutron
        break;

      default:
        E = INF;
    }
  }
  else if ((species == "V") && (Etype == "binding"))
  {
    Real Ef=1.77,gamma=1.0/1.6022e-7;//ev/um^2
    Real r = pow(s*Vatom*3/4/PI,1.0/3); //cluster effective radius
    return Ef-2*gamma*Vatom/r;
  }
  else if ((species == "I") && (Etype == "binding"))
  {
    Real Ef=1.77,gamma=1.0/1.6022e-7;//ev/um^2
    Real r = pow(s*Vatom*3/4/PI,1.0/3); //cluster effective radius
    return Ef-2*gamma*Vatom/r;
  }
  else
    mooseError("Energy not defined for " + Etype + " " + species);
  return E;
}

Real GroupingTest::D_prefactor(int /*s*/, std::string /*species*/) const
{
  Real D0 = 1.0e6;//um^2/s
  return D0;
}

//size S1 and S2
Real GroupingTest::absorb(int S1, int S2, std::string C1, std::string C2,Real T, int tag1, int tag2) const{
    if(tag1==0 && tag2==0) return 0.0;//tag1, tag2 denotes the mobility of C1 and C2; 1: mobile, 0: immobile
    int S = (S1>S2)? S1: S2;
    Real w = pow(48.0*PI*PI/Vatom/Vatom*S,1.0/3);
    Real D_s1 = D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
    Real D_s2 = D_prefactor(S2,C2)*exp(-energy(S2,C2,"migration")/Boltz_const/T);
    return w*Vatom*(D_s1*tag1+D_s2*tag2);//add Vatom for unit concern,  P5/19 in ref

//test the other expression, larger than the one above (3 times)
/*
    Real r_vi = 0.65e-3;//recombination radius in um
    Real r1 = pow(S1*Vatom*3/4/PI,1.0/3); //cluster effective radius
    Real r2 = pow(S2*Vatom*3/4/PI,1.0/3); //cluster effective radius
    return 4*PI*(D_s1*tag1+D_s2*tag2)*(r1+r2+r_vi);
*/
}

Real GroupingTest::diff(int S1, std::string C1,Real T) const {
	return D_prefactor(S1,C1)*exp(-energy(S1,C1,"migration")/Boltz_const/T);
}//in um^2/s

Real
GroupingTest::emit(int S1, int S2, Real T, std::string C1, std::string /*C2*/, int tag1, int tag2) const
{
  //for now only consider self species emission, S1 emits S2, S1==1
  Real emit_c = 0.0;
  if (S1 > S2 && S2==1)
    emit_c = absorb(S1,S2,C1,C1,T,tag1,tag2)/(Vatom) *exp(-energy(S1,C1,"binding")/Boltz_const/T);

  // unit:/s only emit point defect of the same species
  return emit_c;
}

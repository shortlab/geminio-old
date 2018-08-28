/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "GMobile.h"
#include "Conversion.h"

registerMooseObject("GeminioApp", GMobile);

template<>
InputParameters validParams<GMobile>()
{
  InputParameters params = validParams<GBase>();
  return params;
}

GMobile::GMobile(const InputParameters & parameters) : GBase(parameters)
{
  _max_i = (_gc.GroupScheme_i.size() > 0 ? _gc.GroupScheme_i.back() : 0);
  _max_v = (_gc.GroupScheme_v.size() > 0 ? _gc.GroupScheme_v.back() : 0);
}

Real
GMobile::computeQpResidual()
{
  Real res_sum = 0.0;
  int cur_size;//should be positive value
  int ii = _max_mobile_i;
  int vv = _max_mobile_v;
  int max_vi;
  Real conc,conci,concj;

  //printf("return mobile initial: %f %d\n",res_sum,_cur_size);

  if (_cur_size > 0)
  {
    //v type
    cur_size = _cur_size;
    max_vi = std::min(_cur_size + ii, _max_v);


    //vi reaction loss(-)
    for (int i = 1;i <= _max_i; ++i){
      conc = getConcBySize(-i);
      res_sum += conc * _u[_qp] * _gc._absorb(cur_size,-i);
      //printf("vi reaction %d (-): %d %d\n",cur_size,cur_size,-i);
    }

    //vv reaction loss(-)
    for (int i = 1; i <= _max_v-cur_size; ++i)
    {
      // guarantee the largest size doesn't exceed _number_v
      //printf("vv reaction %d (-): %d %d\n",cur_size,cur_size,i);
      conc = getConcBySize(i);
      res_sum += conc*_u[_qp]*_gc._absorb(cur_size,i);
    }
    if(cur_size*2 <= _max_v){
      //printf("vv reaction %d (-): %d %d\n",cur_size,cur_size,cur_size);
      res_sum += _u[_qp]*_u[_qp]*_gc._absorb(cur_size,cur_size);
    }


    //vv reaction gain(+)
    for (int i=1;i <= (int)(cur_size/2);i++)
    {
      //printf("vv reaction %d (+): %d %d\n",cur_size,cur_size-i,cur_size);
      conci = getConcBySize(cur_size-i);
      concj = getConcBySize(i);
      res_sum -= conci * concj *_gc._absorb(cur_size-i,i);
    }

    //vi reaction gain(+)
    for(int i=cur_size+1;i<=max_vi;i++)
    {
      // make sure one is mobile
      if (i-cur_size <= ii || i <= vv )
      {
        conci = getConcBySize(cur_size-i);
        concj = getConcBySize(i);
        res_sum -= conci * concj * _gc._absorb(i,cur_size-i);
        //printf("vi reaction %d (+): %d %d\n",cur_size,cur_size-i,i);
      }
    }

    // v emission loss(-)
    if (cur_size!=1)
      res_sum += _u[_qp]*_gc._emit(cur_size);
      //printf("emission loss %d (-): %d\n",cur_size,cur_size);

    //v+1 emission gain(+)
    if (cur_size<_max_v){
      //printf("emission gain %d (+): %d\n",cur_size,cur_size+1);
      conc = getConcBySize(cur_size+1);
      res_sum -= conc *_gc._emit(cur_size+1);
    }

    if (cur_size == 1)
    {
      for (int i = 2; i <= _max_v; ++i)
      {
        //printf("emission gain %d (+): %d\n",cur_size,i);
        conc = getConcBySize(i);
        res_sum -= conc * _gc._emit(i);
        //printf("res: %f %d %f\n",res_sum,cur_size,_gc._emit(i));
      }
    }

    //dislocation loss(-)
    res_sum += _u[_qp]*_gc._disl(cur_size);
  }

  // i type
  else{
    cur_size = -_cur_size;//make it positive
    max_vi = std::min(cur_size+vv,_max_i);

    //iv reaction loss(-)
    for(int i=1;i<=_max_v;i++){
      conc = getConcBySize(i);
      res_sum += conc *_u[_qp]*_gc._absorb(-cur_size,i);
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,i);
    }

    //ii reaction loss(-)
    for(int i=1;i <= _max_i-cur_size;i++){//guarantee the largest size doesn't exceed _number_v
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,-i);
      conc = getConcBySize(-i);
      res_sum += conc *_u[_qp]*_gc._absorb(-cur_size,-i);
    }
    if(cur_size*2<=_max_i){
      res_sum += _u[_qp]*_u[_qp]*_gc._absorb(-cur_size,-cur_size);
      //printf("reaction %d (-): %d %d\n",_cur_size,_cur_size,_cur_size);
    }

    //ii reaction gain(+)
    for(int i=1;i <= (int)(cur_size/2);i++){
        conci = getConcBySize(-i);
        concj = getConcBySize(i-cur_size);
        res_sum -= conci * concj *_gc._absorb(i-cur_size,-i);
        //printf("reaction %d (+): %d %d\n",_cur_size,i-cur_size,-i);
      //}
    }

    //iv reaction gain(+)
    for(int i=cur_size+1;i<=max_vi;i++){
      if(i-cur_size <= vv || i <= ii){//make sure one is mobile
        conci = getConcBySize(-i);
        concj = getConcBySize(i-cur_size);
        res_sum -= conci * concj *_gc._absorb(-i,i-cur_size);
        //printf("reaction %d (+): %d %d\n",_cur_size,i-cur_size,-i);
      }
    }

    //i emission loss(-)
    if(cur_size!=1){
      res_sum += _u[_qp]*_gc._emit(-cur_size);
      //printf("emission %d (-): %d\n",_cur_size,_cur_size);
    }

    //i+1 emission gain(+)
    if(cur_size<_max_i){
      conc = getConcBySize(-cur_size-1);
      res_sum -= conc *_gc._emit(-cur_size-1);
      //printf("emission %d (+): %d\n",_cur_size,-cur_size-1);
    }
    if(cur_size==1){
      for(int i=2;i<=_max_i;i++){
        conc = getConcBySize(-i);
        res_sum -= conc *_gc._emit(-i);
       // printf("emission %d (+): %d\n",_cur_size,-i);
      }
    }

    //dislocation loss(-)
    res_sum += _u[_qp]*_gc._disl(-cur_size);
  }
  return res_sum*_test[_i][_qp];
}

Real
GMobile::computeQpJacobian()
{
  Real jac_sum = 0.0;
  int cur_size;//should be positive value
  Real conc;
  if(_cur_size>0){//v type

    cur_size = _cur_size;

    //vi reaction loss(-)
    for(int i=1;i<=_max_i;i++){
      conc = getConcBySize(-i);
      jac_sum += conc *_gc._absorb(cur_size,-i);
    }

    //vv reaction loss(-)
    for(int i=1;i <= _max_v-cur_size;i++){//guarantee the largest size doesn't exceed _number_v
      conc = getConcBySize(i);
      jac_sum += conc *_gc._absorb(cur_size,i);//10 is test
    }
    if(cur_size*2<=_max_v)//2*u^2->4*u*phi
      jac_sum += 3.0*_u[_qp]*_gc._absorb(cur_size,cur_size);
    //(*_val_v_vars[cur_size-1])[_qp]

    //v emission loss(-)
    jac_sum += _gc._emit(cur_size);

    //dislocation loss(-)
    jac_sum += _gc._disl(cur_size);
  }

  // i type
  else{
    cur_size = -_cur_size;//make it positive

    //iv reaction loss(-)
    for(int i=1;i<=_max_v;i++){
      conc = getConcBySize(i);
      jac_sum += conc *_gc._absorb(-cur_size,i);
    }

    //ii reaction loss(-)
    for(int i=1;i <= _max_i-cur_size;i++)
    {
      // guarantee the largest size doesn't exceed _number_v
      conc = getConcBySize(-i);
      jac_sum += conc*_gc._absorb(-cur_size,-i);
    }
    if (cur_size * 2 <= _max_i)
      jac_sum += 3.0*_u[_qp]*_gc._absorb(-cur_size,-cur_size);// *(*_val_i_vars[cur_size-1])[_qp]

    //i emission loss(-)
    jac_sum += _gc._emit(-cur_size);

    //dislocation loss(-)
    jac_sum += _gc._disl(-cur_size);
  }
  return jac_sum*_test[_i][_qp] * _phi[_j][_qp];
}

double
GMobile::getConcBySize(int i)
{
  if (i > 0)
  {
    int g = _gc.CurrentGroupV(i);
    return (*_val_v_vars[2*(g-1)])[_qp]+(*_val_v_vars[2*(g-1)+1])[_qp]*(i-_gc.GroupScheme_v_avg[g-1]);
  }
  else
  {
    int g = _gc.CurrentGroupI(-i);
    return (*_val_i_vars[2*(g-1)])[_qp]+(*_val_i_vars[2*(g-1)+1])[_qp]*(-i-_gc.GroupScheme_i_avg[g-1]);
  }
}

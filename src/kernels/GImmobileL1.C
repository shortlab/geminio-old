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

#include "GImmobileL1.h"
#include "Conversion.h"

template<>
InputParameters validParams<GImmobileL1>()
{
  InputParameters params = validParams<GBase>();
  return params;
}

GImmobileL1::GImmobileL1(const InputParameters & parameters) : GBase(parameters)
{
  if (_cur_size>0 && _cur_size< _max_mobile_v)
    mooseError("Please check the mobile vacancy cluster size range, should be immobile>mobile");
  else if (_cur_size<0 && -_cur_size<_max_mobile_i)
    mooseError("Please check the mobile interstitial cluster size range, should be immobile>mobile");
}

Real
GImmobileL1::computeQpResidual()
{
  int cur_size,other_size,group_num;
  Real res_sum = 0.0;
  Real conc1,conc2,conc;
  //printf("return immobile initial: %f %d\n",res_sum,_cur_size);

  // v type
  if (_cur_size > 0)
  {
    //printf("L1: Start; current V group: %d\n",_cur_size);
    cur_size = _cur_size;
    if (_gc.GroupScheme_v_sq[cur_size-1]< 1.0e-12)
      return 0.0;

    int index = 2*_max_mobile_v;//current variable index (start from 0)
    double coefi_1 = (-1-_gc.GroupScheme_v_del[cur_size-1])/2.0;
    double coefi = (-1+_gc.GroupScheme_v_del[cur_size-1])/2.0;

    //left boundary x_{i-1}+1, absorb the same species
    for(int i=0;i<=_max_mobile_v-1;i++){
      int tmp_size = std::min(_max_mobile_v-1-i,_gc.GroupScheme_v_del[cur_size-1]-1);
      group_num =  _gc.CurrentGroupV(_gc.GroupScheme_v[cur_size-1]-i);
      other_size = _max_mobile_v+_max_mobile_v-(cur_size-group_num);
      conc1 = (*_val_v_vars[2*other_size])[_qp]+(*_val_v_vars[2*other_size+1])[_qp]*(_gc.GroupScheme_v[cur_size-1]-i-_gc.GroupScheme_v_avg[group_num-1]);
      for(int j=0;j<=tmp_size;j++){
        conc2 = (*_val_v_vars[2*(i+j)])[_qp];
        res_sum -= (coefi_1+j+1)*conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size-1]-i,i+j+1);
        //printf("absorb %d (vv gain): %d and %d; var: %d %d\n",_cur_size,_gc.GroupScheme_v[cur_size-1]-i,i+j+1,2*(i+j),2*other_size);
      }//vv (gain)
    }

    //left boundary x_{i-1}+1, emission
    conc = (*_val_v_vars[2*index])[_qp]+(_gc.GroupScheme_v[cur_size-1]+1-_gc.GroupScheme_v_avg[cur_size-1])*_u[_qp];
    res_sum += (coefi_1+1)*conc * _gc._emit(_gc.GroupScheme_v[cur_size-1]+1);//v emit (loss)
    //printf("emit %d (v loss): %d; var: %d\n",_cur_size,_gc.GroupScheme_v[cur_size-1]+1,2*index);//v emit (loss)

    //left boundary x_{i-1}+1, absorb the opposite species
    for(int i=0;i<=_max_mobile_i-1;i++){
      int tmp_size = std::min(_max_mobile_i-1-i,_gc.GroupScheme_v_del[cur_size-1]-1);
      for(int j=0;j<=tmp_size;j++){
        conc1 = (*_val_v_vars[2*index])[_qp] + (_gc.GroupScheme_v[cur_size-1]+j+1-_gc.GroupScheme_v_avg[cur_size-1])*_u[_qp];
        conc2 = (*_val_i_vars[2*(i+j)])[_qp];
        res_sum += (coefi_1+j+1)*conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size-1]+j+1,-(i+j+1));
        //printf("absorb %d (vi loss): %d and %d; var: v%d i%d\n",_cur_size,_gc.GroupScheme_v[cur_size-1]+j+1,-(i+j+1),2*index,2*(i+j));
      }//vi (loss)
    }


    if (cur_size != (int)(_gc.GroupScheme_v.size()-1)){

      //right boundary x_{i}+1, absorb the same species
      int tmp_size = std::min(_max_mobile_v-1,_gc.GroupScheme_v_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        for(int j=0;j<=_max_mobile_v-1-i;j++){
          conc1 = (*_val_v_vars[2*index])[_qp]+ (_gc.GroupScheme_v[cur_size]-i-_gc.GroupScheme_v_avg[cur_size-1])*_u[_qp] ;
          conc2 = (*_val_v_vars[2*(i+j)])[_qp];
          res_sum += (coefi-i)*conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size]-i,i+j+1);
          //printf("absorb %d (vv loss): %d and %d; var: %d %d\n",_cur_size,_gc.GroupScheme_v[cur_size]-i,i+j+1,2*(i+j),2*index);
        }//vv (loss)
      }

      //right boundary x_{i}+1, absorb the opposite species
      tmp_size = std::min(_max_mobile_i-1,_gc.GroupScheme_v_del[cur_size-1]-1);
      for (int i=0;i<=tmp_size;i++){
        int tmp2 = std::min(_max_mobile_i-1-i,_gc.GroupScheme_v.back()-_gc.GroupScheme_v[cur_size]-1);
        for (int j=0;j<=tmp2;j++){
          group_num =  _gc.CurrentGroupV(_gc.GroupScheme_v[cur_size]+j+1);
          //other_size = index+i+j+1;
          other_size = index+(group_num-cur_size);
          conc1 = (*_val_v_vars[2*other_size])[_qp]+(*_val_v_vars[2*other_size+1])[_qp]*(_gc.GroupScheme_v[cur_size]+j+1-_gc.GroupScheme_v_avg[group_num-1]);
          conc2 = (*_val_i_vars[2*(i+j)])[_qp];
          res_sum -= (coefi-i)*conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size]+j+1,-(i+j+1));
          //printf("absorb %d (vi gain): %d and %d; var: v%d i%d\n",_cur_size,_gc.GroupScheme_v[cur_size]+j+1,-(i+j+1),2*other_size,2*(i+j));
        }//vi (gain)
      }

      //right boundary x_{i}+1, emission
      conc = (*_val_v_vars[2*index+2])[_qp]+(_gc.GroupScheme_v[cur_size]+1-_gc.GroupScheme_v_avg[cur_size])*(*_val_v_vars[2*index+3])[_qp];
      res_sum -= coefi * conc * _gc._emit(_gc.GroupScheme_v[cur_size]+1);//v emit (gain)
      //printf("emit %d (v gain): %d; var: %d\n",_cur_size,_gc.GroupScheme_v[cur_size]+1,2*(index+1));
    }

    //inside interval
    for(int k=_gc.GroupScheme_v[cur_size-1]+1;k<=_gc.GroupScheme_v[cur_size];k++){
      int tmp_size = std::min(_max_mobile_v,_gc.GroupScheme_v[cur_size]-k);
      conc1 = (*_val_v_vars[2*index])[_qp]+ (k-_gc.GroupScheme_v_avg[cur_size-1])*_u[_qp] ;
      for(int j=1;j<=tmp_size;j++){
        conc2 = (*_val_v_vars[2*(j-1)])[_qp];
        res_sum -= conc1 * conc2 * _gc._absorb(k,j) * j;
        //printf("absorb %d (vv loss): %d and %d; var: %d %d\n",_cur_size,k,j,2*(j-1),2*index);
      }//vv
      tmp_size = std::min(_max_mobile_i,k-_gc.GroupScheme_v[cur_size-1]-1);
      for(int j=1;j<=tmp_size;j++){
        conc2 = (*_val_i_vars[2*(j-1)])[_qp];
        res_sum -= conc1 * conc2 * _gc._absorb(k,-j)*(-j);
        //printf("absorb %d (vi loss): %d and %d; var: v%d i%d\n",_cur_size,k,-j,2*index,2*(j-1));
      }//vi
      res_sum += conc1*_gc._emit(k);//need make up for the beginning point
      //printf("emit %d (v loss): %d; var: %d\n",_cur_size,k,2*index);//v emit (loss)

    }
    conc1 = (*_val_v_vars[2*index])[_qp] + (_gc.GroupScheme_v[cur_size-1]+1-_gc.GroupScheme_v_avg[cur_size-1])*_u[_qp];
    res_sum -= conc1*_gc._emit(_gc.GroupScheme_v[cur_size-1]+1);//makeup
    //printf("emit %d (v gain): %d; var: %d\n",_cur_size,_gc.GroupScheme_v[cur_size-1]+1,2*index);

    //if(res_sum>1.0e-10) printf("return immobile residual L1: %.9f %d\n",res_sum,_cur_size);

    //printf("L1 END\n");
    return 1.0/(_gc.GroupScheme_v_del[cur_size-1]*_gc.GroupScheme_v_sq[cur_size-1])*res_sum *_test[_i][_qp];

  }

  else{
    cur_size = -_cur_size;
    if(_gc.GroupScheme_i_sq[cur_size-1]< 1.0e-12) return 0.0;
    int index = 2*_max_mobile_i;//current variable index (start from 0)
    double coefi_1 = (-1-_gc.GroupScheme_i_del[cur_size-1])/2.0;
    double coefi = (-1+_gc.GroupScheme_i_del[cur_size-1])/2.0;

    //left boundary x_{i-1}+1, absorb the same species
    for(int i=0;i<=_max_mobile_i-1;i++){
      int tmp_size = std::min(_max_mobile_i-1-i,_gc.GroupScheme_i_del[cur_size-1]-1);
      group_num =  _gc.CurrentGroupI(_gc.GroupScheme_i[cur_size-1]-i);
      other_size = _max_mobile_i+_max_mobile_i-(cur_size-group_num);
      conc1 = (*_val_i_vars[2*other_size])[_qp]+(*_val_i_vars[2*other_size+1])[_qp]*(_gc.GroupScheme_i[cur_size-1]-i-_gc.GroupScheme_i_avg[group_num-1]);
      for(int j=0;j<=tmp_size;j++){
        conc2 = (*_val_i_vars[2*(i+j)])[_qp];
        res_sum -= (coefi_1+j+1)*conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size-1]-i),-(i+j+1));
        //printf("ii gain, conc1: %f, conc2: %f, size: %d %d, absorb %f\n",conc1,conc2,-(_gc.GroupScheme_i[cur_size-1]-i),-(i+j+1),_gc._absorb(-(_gc.GroupScheme_i[cur_size-1]-i),-(i+j+1)));
      }//ii (gain)
    }

    //left boundary x_{i-1}+1, emission
    conc = (*_val_i_vars[2*index])[_qp]+(_gc.GroupScheme_i[cur_size-1]+1-_gc.GroupScheme_i_avg[cur_size-1])*_u[_qp];
    res_sum += (coefi_1+1)*conc * _gc._emit(-(_gc.GroupScheme_i[cur_size-1]+1));//i emit (loss)
    //printf("i emit loss, conc: %f, size: %d\n",conc,-(_gc.GroupScheme_i[cur_size-1]+1));//v emit (loss)

    //left boundary x_{i-1}+1, absorb the opposite species
    for(int i=0;i<=_max_mobile_v-1;i++){
      int tmp_size = std::min(_max_mobile_v-1-i,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int j=0;j<=tmp_size;j++){
        conc1 = (*_val_i_vars[2*index])[_qp]+ (_gc.GroupScheme_i[cur_size-1]+j+1-_gc.GroupScheme_i_avg[cur_size-1])*_u[_qp] ;
        conc2 = (*_val_v_vars[2*(i+j)])[_qp];
        res_sum += (coefi_1+j+1)*conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size-1]+j+1),(i+j+1));
      }//iv (loss)
    }


    if(cur_size != (int)(_gc.GroupScheme_i.size()-1)){

      //right boundary x_{i}+1, absorb the same species
      int tmp_size = std::min(_max_mobile_i-1,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        for(int j=0;j<=_max_mobile_i-1-i;j++){
          conc1 = (*_val_i_vars[2*index])[_qp]+ (_gc.GroupScheme_i[cur_size]-i-_gc.GroupScheme_i_avg[cur_size-1])*_u[_qp] ;
          conc2 = (*_val_i_vars[2*(i+j)])[_qp];
          res_sum += (coefi-i)*conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size]-i),-(i+j+1));
        }//ii (loss)
      }

      //right boundary x_{i}+1, absorb the opposite species
      tmp_size = std::min(_max_mobile_v-1,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        int tmp2 = std::min(_max_mobile_v-1-i,_gc.GroupScheme_i.back()-_gc.GroupScheme_i[cur_size]-1);
        for(int j=0;j<=tmp2;j++){
          group_num =  _gc.CurrentGroupI(_gc.GroupScheme_i[cur_size]+j+1);
          //other_size = index+i+j+1;
          other_size = index+(group_num-cur_size);
          conc1 = (*_val_i_vars[2*other_size])[_qp]+(*_val_i_vars[2*other_size+1])[_qp]*(_gc.GroupScheme_i[cur_size]+j+1-_gc.GroupScheme_i_avg[group_num-1]);
          conc2 = (*_val_v_vars[2*(i+j)])[_qp];
          res_sum -= (coefi-i)*conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size]+j+1),(i+j+1));
        }//iv (gain)
      }

      //right boundary x_{i}+1, emission
      conc = (*_val_i_vars[2*index+2])[_qp]+(_gc.GroupScheme_i[cur_size]+1-_gc.GroupScheme_i_avg[cur_size])*(*_val_i_vars[2*index+3])[_qp];
      res_sum -= coefi * conc * _gc._emit(-(_gc.GroupScheme_i[cur_size]+1));//i emit (gain)
    }

    //inside interval
    for(int k=_gc.GroupScheme_i[cur_size-1]+1;k<=_gc.GroupScheme_i[cur_size];k++){
      int tmp_size = std::min(_max_mobile_i,_gc.GroupScheme_i[cur_size]-k);
      conc1 = (*_val_i_vars[2*index])[_qp]+ (k-_gc.GroupScheme_i_avg[cur_size-1])*_u[_qp] ;
      for(int j=1;j<=tmp_size;j++){
        conc2 = (*_val_i_vars[2*(j-1)])[_qp];
        res_sum -= conc1 * conc2 * _gc._absorb(-k,-j) * j;
      }//ii
      tmp_size = std::min(_max_mobile_v,k-_gc.GroupScheme_i[cur_size-1]-1);
      for(int j=1;j<=tmp_size;j++){
        conc2 = (*_val_v_vars[2*(j-1)])[_qp];
        res_sum -= conc1 * conc2 * _gc._absorb(-k,j)*(-j);
      }//iv
      res_sum += conc1*_gc._emit(-k);//need make up for the beginning point
    }
    conc1 = (*_val_i_vars[2*index])[_qp]+ (_gc.GroupScheme_i[cur_size-1]+1-_gc.GroupScheme_i_avg[cur_size-1])*_u[_qp] ;
    res_sum -= conc1*_gc._emit(-(_gc.GroupScheme_i[cur_size-1]+1));//makeup

    //if(res_sum>1.0e-10) printf("return immobile residual L1: %.9f %d\n",res_sum,_cur_size);

    return 1.0/(_gc.GroupScheme_i_del[cur_size-1]*_gc.GroupScheme_i_sq[cur_size-1])*res_sum *_test[_i][_qp];
  }
}

Real
GImmobileL1::computeQpJacobian()
{
  int cur_size;
  Real jac_sum = 0.0;
  Real conc,conc1,conc2;

  if (_cur_size>0){//v type
    cur_size = _cur_size;
    if(_gc.GroupScheme_v_sq[cur_size-1]< 1.0e-12) return 0.0;

    double coefi_1 = (-1-_gc.GroupScheme_v_del[cur_size-1])/2.0;
    double coefi = (-1+_gc.GroupScheme_v_del[cur_size-1])/2.0;

    //left boundary x_{i-1}+1, emission
    conc = _gc.GroupScheme_v[cur_size-1]+1-_gc.GroupScheme_v_avg[cur_size-1];
    jac_sum += (coefi_1+1)*conc * _gc._emit(_gc.GroupScheme_v[cur_size-1]+1);//v emit (loss)

    //left boundary x_{i-1}+1, absorb the opposite species
    for(int i=0;i<=_max_mobile_i-1;i++){
      int tmp_size = std::min(_max_mobile_i-1-i,_gc.GroupScheme_v_del[cur_size-1]-1);
      for(int j=0;j<=tmp_size;j++){
        conc1 = _gc.GroupScheme_v[cur_size-1]+j+1-_gc.GroupScheme_v_avg[cur_size-1];
        conc2 = (*_val_i_vars[2*(i+j)])[_qp];
        jac_sum += (coefi_1+j+1)*conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size-1]+j+1,-(i+j+1));
      }//vi (loss)
    }


    if (cur_size != (int)(_gc.GroupScheme_v.size()-1)){

      //right boundary x_{i}+1, absorb the same species
      int tmp_size = std::min(_max_mobile_v-1,_gc.GroupScheme_v_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        for(int j=0;j<=_max_mobile_v-1-i;j++){
          conc1 = _gc.GroupScheme_v[cur_size]-i-_gc.GroupScheme_v_avg[cur_size-1];
          conc2 = (*_val_v_vars[2*(i+j)])[_qp];
          jac_sum += (coefi-i)*conc1 * conc2 * _gc._absorb(_gc.GroupScheme_v[cur_size]-i,i+j+1);
        }//vv (loss)
      }
    }

    //inside interval
    for(int k=_gc.GroupScheme_v[cur_size-1]+1;k<=_gc.GroupScheme_v[cur_size];k++){
      int tmp_size = std::min(_max_mobile_v,_gc.GroupScheme_v[cur_size]-k);
      conc1 = k-_gc.GroupScheme_v_avg[cur_size-1];
      for(int j=1;j<=tmp_size;j++){
        conc2 = (*_val_v_vars[2*(j-1)])[_qp];
        jac_sum -= conc1 * conc2 * _gc._absorb(k,j) * j;
      }//vv
      tmp_size = std::min(_max_mobile_i,k-_gc.GroupScheme_v[cur_size-1]-1);
      for (int j=1;j<=tmp_size;j++){
        conc2 = (*_val_i_vars[2*(j-1)])[_qp];
        jac_sum -= conc1 * conc2 * _gc._absorb(k,-j)*(-j);
      }//vi
      jac_sum += conc1*_gc._emit(k);//need make up for the beginning point
    }
    conc1 = _gc.GroupScheme_v[cur_size-1]+1-_gc.GroupScheme_v_avg[cur_size-1];
    jac_sum -= conc1*_gc._emit(_gc.GroupScheme_v[cur_size-1]+1);//makeup

    //if(jac_sum>1.0e-10) printf("return immobile gradient L1: %.9f %d\n",jac_sum,_cur_size);
    return 1.0/(_gc.GroupScheme_v_del[cur_size-1]*_gc.GroupScheme_v_sq[cur_size-1])*jac_sum *_test[_i][_qp]*_phi[_j][_qp];
  }

  else{
    cur_size = -_cur_size;
    if(_gc.GroupScheme_i_sq[cur_size-1]< 1.0e-12) return 0.0;

    double coefi_1 = (-1-_gc.GroupScheme_i_del[cur_size-1])/2.0;
    double coefi = (-1+_gc.GroupScheme_i_del[cur_size-1])/2.0;


    //left boundary x_{i-1}+1, emission
    conc = _gc.GroupScheme_i[cur_size-1]+1-_gc.GroupScheme_i_avg[cur_size-1];
    jac_sum += (coefi_1+1)*conc * _gc._emit(-(_gc.GroupScheme_i[cur_size-1]+1));//i emit (loss)

    //left boundary x_{i-1}+1, absorb the opposite species
    for (int i=0;i<=_max_mobile_v-1;i++){
      int tmp_size = std::min(_max_mobile_v-1-i,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int j=0;j<=tmp_size;j++){
        conc1 = _gc.GroupScheme_i[cur_size-1]+j+1-_gc.GroupScheme_i_avg[cur_size-1];
        conc2 = (*_val_v_vars[2*(i+j)])[_qp];
        jac_sum += (coefi_1+j+1)*conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size-1]+j+1),(i+j+1));
      }//iv (loss)
    }

    if (cur_size != (int)(_gc.GroupScheme_i.size()-1)){

      //right boundary x_{i}+1, absorb the same species
      int tmp_size = std::min(_max_mobile_i-1,_gc.GroupScheme_i_del[cur_size-1]-1);
      for(int i=0;i<=tmp_size;i++){
        for(int j=0;j<=_max_mobile_i-1-i;j++){
          conc1 = _gc.GroupScheme_i[cur_size]-i-_gc.GroupScheme_i_avg[cur_size-1];
          conc2 = (*_val_i_vars[2*(i+j)])[_qp];
          jac_sum += (coefi-i)*conc1 * conc2 * _gc._absorb(-(_gc.GroupScheme_i[cur_size]-i),-(i+j+1));
        }//ii (loss)
      }
    }

    //inside interval
    for (int k=_gc.GroupScheme_i[cur_size-1]+1;k<=_gc.GroupScheme_i[cur_size];k++){
      int tmp_size = std::min(_max_mobile_i,_gc.GroupScheme_i[cur_size]-k);
      conc1 = k-_gc.GroupScheme_i_avg[cur_size-1];
      for (int j=1;j<=tmp_size;j++){
        conc2 = (*_val_i_vars[2*(j-1)])[_qp];
        jac_sum -= conc1 * conc2 * _gc._absorb(-k,-j) * j;
      }//ii
      tmp_size = std::min(_max_mobile_v,k-_gc.GroupScheme_i[cur_size-1]-1);
      for (int j=1;j<=tmp_size;j++){
        conc2 = (*_val_v_vars[2*(j-1)])[_qp];
        jac_sum -= conc1 * conc2 * _gc._absorb(-k,j)*(-j);
      }//iv
      jac_sum += conc1*_gc._emit(-k);//need make up for the beginning point
    }
    conc1 = _gc.GroupScheme_i[cur_size-1]+1-_gc.GroupScheme_i_avg[cur_size-1];
    jac_sum -= conc1*_gc._emit(-(_gc.GroupScheme_i[cur_size-1]+1));//makeup


    //if(jac_sum>1.0e-10) printf("return immobile gradient L1: %.9f %d\n",jac_sum,_cur_size);
    return 1.0/(_gc.GroupScheme_i_del[cur_size-1]*_gc.GroupScheme_i_sq[cur_size-1])*jac_sum *_test[_i][_qp]*_phi[_j][_qp];
  }
}

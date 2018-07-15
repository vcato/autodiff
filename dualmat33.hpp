#ifndef DUALMAT33_HPP
#define DUALMAT33_HPP

#include "mat33.hpp"
#include "dualfloat.hpp"

namespace autodiff {


using DualMat33 = Mat33<DualFloat>;


inline DualMat33 dual(Mat33<float> v,Mat33<float> &dv)
{
  DualFloat values[3][3] = {
    {dual(v[0][0],dv[0][0]),dual(v[0][1],dv[0][1]),dual(v[0][2],dv[0][2])},
    {dual(v[1][0],dv[1][0]),dual(v[1][1],dv[1][1]),dual(v[1][2],dv[1][2])},
    {dual(v[2][0],dv[2][0]),dual(v[2][1],dv[2][1]),dual(v[2][2],dv[2][2])},
  };

  return DualMat33(values);
}


}


#endif // DUALMAT33_HPP

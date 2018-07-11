#include "dualfloat.hpp"
using DualVec3 = Vec3<DualFloat>;


inline DualVec3 dual(FloatVec3 v,FloatVec3 &dv)
{
  return {{v.x(),dv.x()},{v.y(),dv.y()},{v.z(),dv.z()}};
}
